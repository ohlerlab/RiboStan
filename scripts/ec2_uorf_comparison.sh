#!/usr/bin/env bash
# PI research EC2: Ribostan vs RiboTaper uORF comparison on chr22.
set -euxo pipefail
export MAMBA_ROOT_PREFIX=/opt/micromamba
eval "$(micromamba shell hook -s bash)"

WORKDIR=/home/ec2-user/work/uorf_compare
RIBOSTAN=/home/ec2-user/work/RiboStan
mkdir -p "$WORKDIR"
cd "$WORKDIR"
exec > "$WORKDIR/run.log" 2>&1

echo "=== $(date -Is) starting uORF comparison ==="

# Bootstrap ran as root; ec2-user needs write access to micromamba cache.
sudo chown -R ec2-user:ec2-user /opt/micromamba 2>/dev/null || true

# Separate env from ribotaper (python=2.7 pin breaks modern Bioconductor solver).
if ! micromamba env list | grep -q '/envs/ribostan'; then
  micromamba create -y -n ribostan -c conda-forge -c bioconda \
    r-base=4.4 r-rcpp r-dplyr r-magrittr r-stringr r-readr r-tidyr r-broom \
    r-tibble r-matrix r-ggplot2 r-rstan r-multitaper r-testthat \
    bioconductor-rtracklayer bioconductor-genomicfeatures bioconductor-genomicranges \
    bioconductor-iranges bioconductor-s4vectors bioconductor-biocgenerics \
    bioconductor-biostrings bioconductor-rsamtools bioconductor-genomeinfodb \
    bioconductor-txdbmaker bioconductor-annotationhub bioconductor-genomicalignments
fi

micromamba activate ribostan
cd "$RIBOSTAN"
R CMD INSTALL --no-multiarch --with-keep.source .

cd "$WORKDIR"
Rscript "$RIBOSTAN/scripts/compare_uorfs_prep.R" "$WORKDIR"

# RiboTaper stays in its own env (python 2.7 + old bedtools).
eval "$(micromamba shell hook -s bash)"
micromamba activate ribotaper

READLENGTHS=$(paste -sd, ribo_readlengths.txt)
OFFSETS=$(paste -sd, ribo_offsets.txt)
ANNOT_DIR="$WORKDIR/ribotaper_annot"
mkdir -p "$ANNOT_DIR"

create_annotations_files.bash "$WORKDIR/test.gc32.gtf" "$WORKDIR/chr22.fa" false false "$ANNOT_DIR"

cp "$WORKDIR/nchr22_genomic.bam" "$WORKDIR/ribo.bam"
cp "$WORKDIR/nchr22_genomic.bam" "$WORKDIR/rna.bam"

cd "$WORKDIR"
Ribotaper.sh ribo.bam rna.bam "$ANNOT_DIR" "$READLENGTHS" "$OFFSETS" 2

echo "=== $(date -Is) RiboTaper finished ==="

micromamba activate ribostan
Rscript "$RIBOSTAN/scripts/compare_uorfs_report.R" "$WORKDIR"
echo "=== $(date -Is) COMPARISON_DONE ==="
