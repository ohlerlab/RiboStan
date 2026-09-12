#!/usr/bin/env bash
# Build RiboTaper-exact genomic P_sites_all BED from the matched primary BAM
# using published HEK293 offsets (Calviello et al. 2016 / ORFcalling README):
#   lengths 26,27,28,29  offsets 12,12,12,12
# Then (optionally) compare multitaper p on those tracks vs Ribostan psites
# for C1 RiboTaper uORFs.
set -euo pipefail

WORKDIR="${1:-.}"
THREADS="${2:-8}"
cd "$WORKDIR"

CONDA_ROOT="${CONDA_ROOT:-/home/ec2-user/miniforge3}"
if [[ -f "$CONDA_ROOT/etc/profile.d/conda.sh" ]]; then
  # shellcheck source=/dev/null
  source "$CONDA_ROOT/etc/profile.d/conda.sh"
  conda activate rnaseq
fi

BAM="${HEK293_GENOME_BAM:-hek293_genome_rrnafilt_primary.bam}"
OUTDIR="${HEK293_RIBOTAPER_PSITES_DIR:-ribotaper_psites}"
mkdir -p "$OUTDIR"
cd "$OUTDIR"

# RiboTaper expects RIBO_best.bam in CWD (primary alignments; we already have -F 0x100)
if [[ ! -f RIBO_best.bam ]]; then
  ln -sf "../$BAM" RIBO_best.bam
fi
# Dummy RNA bam so the stock script can finish; we only need P_sites_all.
if [[ ! -f RNA_best.bam ]]; then
  # Create a tiny empty BAM with same header for RNA step (script always runs it).
  samtools view -H RIBO_best.bam | samtools view -b -o RNA_best.bam -
fi

BEDTOOLS_DIR="$(dirname "$(command -v bedtools)")"
LENGTHS="${HEK293_RT_LENGTHS:-26,27,28,29}"
OFFSETS="${HEK293_RT_OFFSETS:-12,12,12,12}"

if [[ ! -f P_sites_all ]]; then
  echo "=== $(date -Is) Computing RiboTaper P_sites_all (lengths=$LENGTHS offsets=$OFFSETS) ==="
  # Inline the P-site half of P_sites_RNA_sites_calc.bash (skip RNA to save time)
  IFS=',' read -r -a lengths <<< "$LENGTHS"
  IFS=',' read -r -a offsets <<< "$OFFSETS"
  if [[ ${#lengths[@]} -ne ${#offsets[@]} ]]; then
    echo "lengths/offsets length mismatch" >&2
    exit 1
  fi
  rm -f tmp_P_sites_* P_sites_all
  for i in "${!lengths[@]}"; do
    len="${lengths[$i]}"
    offs="${offsets[$i]}"
    echo "processing ${len} nt reads with offset +${offs}"
    "$BEDTOOLS_DIR/bamToBed" -cigar -bed12 -i RIBO_best.bam \
      | awk -v env_var="$len" -F'\t' '{split($11,c,","); if((c[1]+c[2]+c[3])==env_var) print $0 "\t" c[1]+c[2]+c[3] > "tmp_align_len"}'
    # shellcheck disable=SC2016
    awk -F'\t' '{split($11,c,","); print $0 "\t" c[1]+c[2]+c[3]+c[4]+c[5] }' tmp_align_len \
      | awk -v env_var="$offs" '{split($11,c,","); split($12,d,","); if($6=="+" && c[1]<env_var && (c[1]+c[2])>env_var) $2=($2)+d[2]+(env_var-c[1]); if($6=="+" && c[1]>=env_var) $2=($2)+env_var; if($6=="+" && c[1]<env_var && (c[1]+c[2])<12 && (c[1]+c[2]+c[3])>env_var) $2=($2)+d[3]+(env_var-c[1]-c[2]); if($6=="+" && (c[1]+c[2]+c[3])<env_var) $2=($2)+d[4]+(env_var-c[1]-c[2]-c[3]); if($6=="-" && c[1]>=$NF-env_var) $2=$2+$NF-env_var-1; if($6=="-" && c[1]<($NF-env_var) && (c[1]+c[2])>=($NF-env_var)) $2=($2)+d[2]+(($NF-env_var-1)-c[1]); if($6=="-" && c[1]<($NF-env_var) && (c[1]+c[2])<($NF-env_var) && (c[1]+c[2]+c[3])>=($NF-env_var)) $2=($2)+d[3]+(($NF-env_var-1)-c[1]-c[2]); if($6=="+" && (c[1]+c[2]+c[3])<($NF-env_var)) $2=($2)+d[4]+(($NF-env_var)-c[1]-c[2]-c[3]); print $0}' OFS='\t' \
      | awk '{$3=$2+1 ; print $0 }' OFS='\t' \
      | awk '{if($2>0) print $0}' OFS='\t' > "tmp_P_sites_${len}"
    rm -f tmp_align_len
  done
  cat tmp_P_sites_* > P_sites_all
  rm -f tmp_P_sites_*
fi

wc -l P_sites_all
echo "Wrote $(pwd)/P_sites_all"
# BED6 subset for faster downstream use
awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$4,$5,$6}' P_sites_all > P_sites_all.bed6
wc -l P_sites_all.bed6
echo "=== $(date -Is) DONE ==="
