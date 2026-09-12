#!/usr/bin/env bash
# Download HEK293 RiboTaper benchmark inputs (Calviello et al. 2016, GSE73136).
set -euxo pipefail

WORKDIR="${1:-.}"
mkdir -p "$WORKDIR"
cd "$WORKDIR"

GENCODE19_BASE="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19"

if [[ ! -f gencode.v19.annotation.gtf.gz ]]; then
  curl -L -o gencode.v19.annotation.gtf.gz \
    "${GENCODE19_BASE}/gencode.v19.annotation.gtf.gz"
fi
if [[ ! -f gencode.v19.transcripts.fa.gz ]] || [[ $(wc -c < gencode.v19.transcripts.fa.gz) -lt 1000000 ]]; then
  rm -f gencode.v19.transcripts.fa.gz gencode.v19.transcripts.fa
  # Release 19 has pc_transcripts + lncRNA, not a single transcripts.fa.gz.
  curl -L -o gencode.v19.transcripts.fa.gz \
    "${GENCODE19_BASE}/gencode.v19.pc_transcripts.fa.gz"
fi
if [[ $(wc -c < gencode.v19.transcripts.fa.gz) -lt 1000000 ]]; then
  echo "ERROR: transcript FASTA download looks truncated ($(wc -c < gencode.v19.transcripts.fa.gz) bytes)" >&2
  exit 1
fi
if [[ ! -f GRCh37.p13.genome.fa.gz ]]; then
  curl -L -o GRCh37.p13.genome.fa.gz \
    "${GENCODE19_BASE}/GRCh37.p13.genome.fa.gz"
fi

gunzip -kf gencode.v19.annotation.gtf.gz
gunzip -kf gencode.v19.transcripts.fa.gz
gunzip -kf GRCh37.p13.genome.fa.gz

if [[ ! -f GSM1887643_HEK293_filtered_ORFs.bed.gz ]]; then
  curl -L -o GSE73136_RAW.tar \
    "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE73136&format=file"
  tar -xf GSE73136_RAW.tar GSM1887643_HEK293_filtered_ORFs.bed.gz
fi

if [[ ! -f HEK293_RIBO.fastq.gz ]] || [[ $(wc -c < HEK293_RIBO.fastq.gz) -lt 100000000 ]]; then
  rm -f HEK293_RIBO.fastq.gz
  if command -v prefetch >/dev/null 2>&1; then
    prefetch SRR2433794
    fasterq-dump --gzip SRR2433794
    mv SRR2433794.fastq.gz HEK293_RIBO.fastq.gz
  else
    # ENA FTP (single-end; no _1 suffix for this run).
    curl -L -o HEK293_RIBO.fastq.gz \
      "https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR243/004/SRR2433794/SRR2433794.fastq.gz"
  fi
fi
if [[ $(wc -c < HEK293_RIBO.fastq.gz) -lt 100000000 ]]; then
  echo "ERROR: Ribo-seq FASTQ download looks truncated ($(wc -c < HEK293_RIBO.fastq.gz) bytes)" >&2
  exit 1
fi

echo "Data ready in $(pwd)"

BENCH="$(cd "$(dirname "$0")" && pwd)"
if [[ ! -f gencode.v19.ccds_appris.gtf ]] || [[ "${HEK293_REBUILD_GTF:-FALSE}" == TRUE ]]; then
  if command -v Rscript >/dev/null 2>&1; then
    echo "Building RiboTaper-style filtered GTF (CCDS+APPRIS) ..."
    Rscript "$BENCH/filter_gtf_ribotaper_style.R" "$WORKDIR"
  else
    echo "NOTE: Rscript not on PATH; run filter_gtf_ribotaper_style.R after 'conda activate ribostan'"
  fi
fi
