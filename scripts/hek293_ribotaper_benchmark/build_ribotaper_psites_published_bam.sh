#!/usr/bin/env bash
# Build RiboTaper P_sites_all from the published ORFcalling rpf_alignment.bam
set -euo pipefail
WORKDIR="${1:-/home/ec2-user/work/hek293_benchmark}"
cd "$WORKDIR"
source /home/ec2-user/miniforge3/etc/profile.d/conda.sh
conda activate rnaseq

PUB_BAM="ribotaper_published_tracks/orfcalling_drive/mapping_files/rpf_alignment.bam"
for i in $(seq 1 40); do
  if [[ -f "$PUB_BAM" ]] && [[ $(stat -c%s "$PUB_BAM") -gt 1000000 ]]; then
    break
  fi
  echo "waiting for $PUB_BAM ..."
  sleep 30
done
ls -lh "$PUB_BAM"

OUT=ribotaper_psites_published_bam
mkdir -p "$OUT"
cd "$OUT"
if [[ ! -f RIBO_best.bam ]]; then
  samtools view -b -F 0x100 -@ 8 "../$PUB_BAM" > RIBO_best.bam
  samtools index RIBO_best.bam
fi
samtools view -H RIBO_best.bam | samtools view -b -o RNA_best.bam -

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
echo PUB_PSITES_DONE
