#!/usr/bin/env bash
# EC2 (miniforge3): Ribostan vs published RiboTaper HEK293 uORFs (GSE73136).
set -euxo pipefail

CONDA_ROOT="${CONDA_ROOT:-/home/ec2-user/miniforge3}"
# shellcheck source=/dev/null
source "$CONDA_ROOT/etc/profile.d/conda.sh"

WORKDIR="${WORKDIR:-/home/ec2-user/work/hek293_benchmark}"
RIBOSTAN="${RIBOSTAN:-/home/ec2-user/work/RiboStan}"
BENCH="${RIBOSTAN}/scripts/hek293_ribotaper_benchmark"
N_CORES="${N_CORES:-$(nproc)}"

mkdir -p "$WORKDIR"
cd "$WORKDIR"
exec > "$WORKDIR/run.log" 2>&1

echo "=== $(date -Is) HEK293 RiboTaper benchmark (conda) ==="

if ! conda env list | grep -q '^ribostan '; then
  conda create -y -n ribostan -c conda-forge -c bioconda \
    r-base=4.4 r-rcpp r-dplyr r-magrittr r-stringr r-readr r-tidyr r-broom \
    r-tibble r-matrix r-ggplot2 r-multitaper r-testthat r-rstan \
    bioconductor-rtracklayer bioconductor-genomicfeatures bioconductor-genomicranges \
    bioconductor-iranges bioconductor-s4vectors bioconductor-biocgenerics \
    bioconductor-biostrings bioconductor-rsamtools bioconductor-genomeinfodb \
    bioconductor-txdbmaker bioconductor-annotationhub
fi

bash "$BENCH/download_data.sh" "$WORKDIR"

conda activate ribostan
cd "$RIBOSTAN"
R CMD INSTALL --no-multiarch --with-keep.source .

cd "$WORKDIR"
export HEK293_N_CORES="$N_CORES"
bash "$BENCH/run_ribotaper_matched.sh" "$WORKDIR" "$N_CORES"

echo "=== $(date -Is) HEK293_COMPARISON_DONE ==="
cat hek293_comparison_summaryribotaper_matched.txt
