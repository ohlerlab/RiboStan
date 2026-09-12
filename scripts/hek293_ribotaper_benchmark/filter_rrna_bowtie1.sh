#!/usr/bin/env bash
# Reproduce RiboTaper's documented rRNA pre-filter for the HEK293 dataset.
# GSM1887643's "Data processing" notes (GEO, GSE73136) state:
#   "Reads aligning to rRNA sequences were discarded using Bowtie (v1.00),
#    with the following parameters: --best. Unaligned Ribo-seq reads and
#    RNA-seq reads were aligned to the human genome..."
# i.e. rRNA removal happens *before* genome alignment, using Bowtie1 --best
# against an rRNA reference, keeping only the reads that do NOT align
# (--un). This step was missing from the first genome-alignment experiment
# (align_genome_ribotaper_style.sh applied directly to the raw FASTQ).
#
# rRNA reference: human 45S precursor (NR_046235.3, covers 18S/5.8S/28S +
# ITS/ETS), 5S (NR_023363.1), and mitochondrial 12S/16S (NC_012920.1
# sub-intervals) - a standard, complete human rRNA target set. RiboTaper's
# own notes do not name the exact rRNA reference/accessions they used, so
# this is the best-available reconstruction, not a byte-identical
# reproduction.
set -euxo pipefail

WORKDIR="${1:-.}"
THREADS="${2:-8}"
RRNA_FASTA="${3:-human_rrna_reference.fa}"
cd "$WORKDIR"

if [[ ! -f "$RRNA_FASTA" ]]; then
  echo "!!! $RRNA_FASTA not found" >&2
  exit 1
fi

if [[ ! -f rrna_index.1.ebwt ]]; then
  bowtie-build "$RRNA_FASTA" rrna_index
fi

bowtie --best -p "$THREADS" -S rrna_index HEK293_RIBO.fastq.gz \
  --un HEK293_RIBO.norrna.fastq /dev/null > bowtie_rrna_run.log 2>&1

echo "=== bowtie rRNA filter summary ==="
tail -5 bowtie_rrna_run.log
echo "=== non-rRNA read count (fastq records) ==="
awk 'END{print NR/4}' HEK293_RIBO.norrna.fastq
