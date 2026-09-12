#!/usr/bin/env bash
# EC2: Ribostan vs published RiboTaper HEK293 uORFs (GSE73136).
set -euxo pipefail
export MAMBA_ROOT_PREFIX=/opt/micromamba
eval "$(micromamba shell hook -s bash)"

WORKDIR=/home/ec2-user/work/hek293_benchmark
RIBOSTAN="${RIBOSTAN:-/home/ec2-user/work/RiboStan}"
BENCH="${RIBOSTAN}/scripts/hek293_ribotaper_benchmark"
mkdir -p "$WORKDIR"
cd "$WORKDIR"
exec > "$WORKDIR/run.log" 2>&1

echo "=== $(date -Is) HEK293 RiboTaper benchmark ==="

sudo chown -R ec2-user:ec2-user /opt/micromamba 2>/dev/null || true

if ! micromamba env list | grep -q '/envs/ribostan'; then
  micromamba create -y -n ribostan -c conda-forge -c bioconda \
    r-base=4.4 r-rcpp r-dplyr r-magrittr r-stringr r-readr r-tidyr r-broom \
    r-tibble r-matrix r-ggplot2 r-rstan r-multitaper r-testthat \
    bioconductor-rtracklayer bioconductor-genomicfeatures bioconductor-genomicranges \
    bioconductor-iranges bioconductor-s4vectors bioconductor-biocgenerics \
    bioconductor-biostrings bioconductor-rsamtools bioconductor-genomeinfodb \
    bioconductor-txdbmaker bioconductor-annotationhub star samtools
fi

micromamba activate ribostan
bash "$BENCH/download_data.sh" "$WORKDIR"
bash "$BENCH/align_transcriptome.sh" "$WORKDIR" "$(nproc)"

cd "$RIBOSTAN"
R CMD INSTALL --no-multiarch --with-keep.source .

cd "$WORKDIR"
Rscript "$BENCH/compare_hek293_prep.R" "$WORKDIR"
Rscript "$BENCH/compare_hek293_report.R" "$WORKDIR"

echo "=== $(date -Is) HEK293_COMPARISON_DONE ==="
cat hek293_comparison_summary.txt
