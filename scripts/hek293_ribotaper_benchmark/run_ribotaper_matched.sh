#!/usr/bin/env bash
# End-to-end Ribostan vs RiboTaper HEK293 comparison on *matched* inputs:
# rRNA pre-filter, genome STAR (RiboTaper filters), primary-only, project to
# transcript space (strand fix), metagene-peak offsets, RiboTaper nested rule,
# raw p<0.05, in-frame gate 0.5. Skips steps when outputs already exist.
set -euo pipefail

WORKDIR="${1:-.}"
THREADS="${2:-$(nproc)}"
RIBOSTAN="${RIBOSTAN:-$(cd "$(dirname "$0")/../.." && pwd)}"
BENCH="${RIBOSTAN}/scripts/hek293_ribotaper_benchmark"
cd "$WORKDIR"

CONDA_ROOT="${CONDA_ROOT:-/home/ec2-user/miniforge3}"
# shellcheck source=/dev/null
source "$CONDA_ROOT/etc/profile.d/conda.sh"

export HEK293_GTF_SOURCE="${HEK293_GTF_SOURCE:-gencode.v19.annotation.gtf}"
export HEK293_USE_CCDS="${HEK293_USE_CCDS:-TRUE}"
export HEK293_USE_APPRIS="${HEK293_USE_APPRIS:-TRUE}"
export HEK293_REBUILD_GTF="${HEK293_REBUILD_GTF:-FALSE}"
# Optional: set HEK293_GTF_NAME to skip auto-build and use a fixed GTF filename.
export HEK293_FDR_CORRECT="${HEK293_FDR_CORRECT:-FALSE}"
export HEK293_NESTED_RULE="${HEK293_NESTED_RULE:-ribotaper}"
export HEK293_MULTIMAP="${HEK293_MULTIMAP:-all}"
export HEK293_MIN_FRAME_FRAC="${HEK293_MIN_FRAME_FRAC:-0.5}"
export HEK293_OFFSET_METHOD="${HEK293_OFFSET_METHOD:-metagene_peak}"
export HEK293_MULTIMAP_QC="${HEK293_MULTIMAP_QC:-FALSE}"
export HEK293_APPLY_PHASE_SHIFT="${HEK293_APPLY_PHASE_SHIFT:-FALSE}"
export HEK293_N_CORES="${HEK293_N_CORES:-$THREADS}"
export HEK293_READ_SOURCE="${HEK293_READ_SOURCE:-genome_primary}"
export HEK293_GENOME_BAM="${HEK293_GENOME_BAM:-hek293_genome_rrnafilt_primary.bam}"
export HEK293_PROJECTED_RDS="${HEK293_PROJECTED_RDS:-rpfs_genome_projected_rrnafilt_matched.rds}"
export HEK293_OUT_SUFFIX="${HEK293_OUT_SUFFIX:-ribotaper_matched}"

GENOME_BAM="$WORKDIR/$HEK293_GENOME_BAM"
PROJECTED_RDS="$WORKDIR/$HEK293_PROJECTED_RDS"

echo "=== $(date -Is) RiboTaper-matched HEK293 benchmark ==="
echo "WORKDIR=$WORKDIR OUT_SUFFIX=$HEK293_OUT_SUFFIX PROJECTED_RDS=$HEK293_PROJECTED_RDS"

if [[ ! -f "$GENOME_BAM" ]]; then
  echo "Genome BAM missing; running rRNA filter + genome alignment ..."
  conda activate rnaseq
  bash "$BENCH/filter_rrna_bowtie1.sh" "$WORKDIR" "$THREADS"
  bash "$BENCH/align_genome_ribotaper_style_rrnafilt.sh" "$WORKDIR" "$THREADS"
fi

if [[ ! -f "$PROJECTED_RDS" ]]; then
  echo "Ensuring RiboTaper-style annotation GTF (CCDS+APPRIS on GENCODE v19) ..."
  Rscript "$BENCH/filter_gtf_ribotaper_style.R" "$WORKDIR"
  echo "Projected RDS missing; running convert_genome_to_txspace.R ..."
  conda activate ribostan
  ulimit -s unlimited
  Rscript "$BENCH/convert_genome_to_txspace.R" "$WORKDIR"
fi

conda activate ribostan
ulimit -s unlimited
Rscript "$BENCH/compare_hek293_prep.R" "$WORKDIR"
Rscript "$BENCH/compare_hek293_report.R" "$WORKDIR"
Rscript "$BENCH/diagnose_differences.R" "$WORKDIR" "$HEK293_OUT_SUFFIX"
Rscript "$BENCH/compare_periodicity_vectors.R" "$WORKDIR" "$HEK293_OUT_SUFFIX"

echo "=== $(date -Is) DONE ==="
cat "$WORKDIR/hek293_comparison_summary${HEK293_OUT_SUFFIX}.txt"
