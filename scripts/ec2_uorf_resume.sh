#!/usr/bin/env bash
# Resume uORF comparison from prep (ribostan env + package already installed).
set -euxo pipefail
export MAMBA_ROOT_PREFIX=/opt/micromamba
eval "$(micromamba shell hook -s bash)"

WORKDIR=/home/ec2-user/work/uorf_compare
RIBOSTAN=/home/ec2-user/work/RiboStan
mkdir -p "$WORKDIR"
cd "$WORKDIR"
exec >> "$WORKDIR/run.log" 2>&1

echo "=== $(date -Is) resuming from prep ==="

ulimit -s unlimited || true

micromamba activate ribostan
if [[ ! -f "$WORKDIR/ribostan_uorfs.rds" || ! -f "$WORKDIR/nchr22_genomic.bam" || ! -f "$WORKDIR/ribo_readlengths.txt" ]]; then
  Rscript "$RIBOSTAN/scripts/compare_uorfs_prep.R" "$WORKDIR"
else
  echo "prep outputs already present; skipping compare_uorfs_prep.R"
fi

eval "$(micromamba shell hook -s bash)"
micromamba activate ribotaper

READLENGTHS=$(paste -sd, ribo_readlengths.txt)
OFFSETS=$(paste -sd, ribo_offsets.txt)
ANNOT_DIR="$WORKDIR/ribotaper_annot"
RIBOTAPER_GTF="$WORKDIR/gencode.v32.chr22.gtf"
mkdir -p "$ANNOT_DIR"

if [[ ! -s "$RIBOTAPER_GTF" ]]; then
  echo "Downloading GENCODE v32 GTF (chr22) for RiboTaper..."
  curl -fsSL "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_32/gencode.v32.annotation.gtf.gz" \
    | gzip -dc | awk '$1 == "chr22"' > "$RIBOTAPER_GTF"
fi

# RiboTaper needs libgfortran.so.3 from the old bioconda stack.
micromamba install -y -n ribotaper -c conda-forge libgfortran-ng || true

create_annotations_files.bash "$RIBOTAPER_GTF" "$WORKDIR/chr22.fa" false false "$ANNOT_DIR"

cp "$WORKDIR/nchr22_genomic.bam" "$WORKDIR/ribo.bam"
cp "$WORKDIR/nchr22_genomic.bam" "$WORKDIR/rna.bam"

cd "$WORKDIR"
Ribotaper.sh ribo.bam rna.bam "$ANNOT_DIR" "$READLENGTHS" "$OFFSETS" 2

echo "=== $(date -Is) RiboTaper finished ==="

micromamba activate ribostan
Rscript "$RIBOSTAN/scripts/compare_uorfs_report.R" "$WORKDIR"
echo "=== $(date -Is) COMPARISON_DONE ==="
