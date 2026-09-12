#!/usr/bin/env bash
# Build RiboTaper-style transcript-space P_sites_all BED from a transcriptome BAM.
# Chrom field = transcript ID (ENST...); coords = 0-based position on transcript.
# Offsets: lengths 26-29 nt, P-site offset +12 (Calviello / ORFcalling README).
set -euo pipefail

WORKDIR="${1:-.}"
OUT_TAG="${2:-published_bam}"
cd "$WORKDIR"

CONDA_ROOT="${CONDA_ROOT:-/home/ec2-user/miniforge3}"
# shellcheck source=/dev/null
source "$CONDA_ROOT/etc/profile.d/conda.sh"
conda activate rnaseq

case "$OUT_TAG" in
  published_bam)
    BAM="ribotaper_published_tracks/orfcalling_drive/mapping_files/rpf_transcript_alignment.bam"
    OUTDIR="ribotaper_psites_transcript_published_bam"
    ;;
  matched_bam)
    BAM="${HEK293_TX_BAM:-hek293_transcript_rrnafilt_primary.bam}"
    OUTDIR="ribotaper_psites_transcript_matched"
    ;;
  *)
    echo "Unknown OUT_TAG: $OUT_TAG (use published_bam or matched_bam)" >&2
    exit 1
    ;;
esac

mkdir -p "$OUTDIR"
cd "$OUTDIR"

if [[ ! -f RIBO_best.bam ]]; then
  ln -sf "../$BAM" RIBO_best.bam
fi
# Transcriptome BAMs often fail samtools index (discontig refs); bamToBed scans
# the full file and does not require an index.

BEDTOOLS_DIR="$(dirname "$(command -v bedtools)")"
rm -f tmp_P_sites_* tmp_align_len P_sites_all P_sites_all.bed6

for pair in 26:12 27:12 28:12 29:12; do
  len="${pair%%:*}"
  offs="${pair##*:}"
  echo "processing ${len} nt offset +${offs}"
  "$BEDTOOLS_DIR/bamToBed" -cigar -bed12 -i RIBO_best.bam \
    | awk -v env_var="$len" -F'\t' '{split($11,c,","); if((c[1]+c[2]+c[3])==env_var) print $0 "\t" c[1]+c[2]+c[3] > "tmp_align_len"}'
  awk -F'\t' '{split($11,c,","); print $0 "\t" c[1]+c[2]+c[3]+c[4]+c[5] }' tmp_align_len \
    | awk -v env_var="$offs" '{split($11,c,","); split($12,d,","); if($6=="+" && c[1]<env_var && (c[1]+c[2])>env_var) $2=($2)+d[2]+(env_var-c[1]); if($6=="+" && c[1]>=env_var) $2=($2)+env_var; if($6=="+" && c[1]<env_var && (c[1]+c[2])<12 && (c[1]+c[2]+c[3])>env_var) $2=($2)+d[3]+(env_var-c[1]-c[2]); if($6=="+" && (c[1]+c[2]+c[3])<env_var) $2=($2)+d[4]+(env_var-c[1]-c[2]-c[3]); if($6=="-" && c[1]>=$NF-env_var) $2=$2+$NF-env_var-1; if($6=="-" && c[1]<($NF-env_var) && (c[1]+c[2])>=($NF-env_var)) $2=($2)+d[2]+(($NF-env_var-1)-c[1]); if($6=="-" && c[1]<($NF-env_var) && (c[1]+c[2])<($NF-env_var) && (c[1]+c[2]+c[3])>=($NF-env_var)) $2=($2)+d[3]+(($NF-env_var-1)-c[1]-c[2]); if($6=="+" && (c[1]+c[2]+c[3])<($NF-env_var)) $2=($2)+d[4]+(($NF-env_var)-c[1]-c[2]-c[3]); print $0}' OFS='\t' \
    | awk '{$3=$2+1 ; print $0 }' OFS='\t' \
    | awk '{if($2>0) print $0}' OFS='\t' > "tmp_P_sites_${len}"
  rm -f tmp_align_len
done

cat tmp_P_sites_* > P_sites_all
awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$4,$5,$6}' P_sites_all > P_sites_all.bed6
rm -f tmp_P_sites_*
wc -l P_sites_all.bed6
echo "TX_PSITES_DONE $OUTDIR"
