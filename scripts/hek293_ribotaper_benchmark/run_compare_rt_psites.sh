#!/usr/bin/env bash
set -euo pipefail
cd /home/ec2-user/work/hek293_benchmark
RSCRIPT=/home/ec2-user/miniforge3/envs/ribostan/bin/Rscript
BEDTOOLS=/home/ec2-user/miniforge3/envs/rnaseq/bin/bedtools
BENCH=/home/ec2-user/work/RiboStan/scripts/hek293_ribotaper_benchmark
ln -sfn "$BEDTOOLS" /tmp/bedtools
export PATH="/home/ec2-user/miniforge3/envs/ribostan/bin:/tmp:$PATH"
ulimit -s unlimited

clean_genome_bed() {
  local src=$1 dest=$2 pre=${dest%.clean.bed6}.prefilter.bed6
  echo "Cleaning genome BED $src"
  rm -f "$pre" "$dest"
  awk 'BEGIN{FS=OFS="\t"} NF==6 && $1 ~ /^chr/ {print $1,$2,$3,".",0,$6}' "$src" > "$pre"
  wc -l "$src" "$pre"
  LC_ALL=C sort -k1,1 -k2,2n --parallel=8 -S 4G "$pre" -o "$dest"
  wc -l "$dest"
  rm -f "$pre"
}

# 1) Build transcript-space P-sites from published ORFcalling BAM (exact inputs)
if [[ ! -f ribotaper_psites_transcript_published_bam/P_sites_all.bed6 ]]; then
  bash "$BENCH/build_ribotaper_psites_transcript.sh" . published_bam
fi

# 2) Transcript-space compare (RiboTaper uORF testing coordinates)
echo "=== Compare: published transcript P-sites ==="
export HEK293_RT_PSITE_MODE=transcript
export HEK293_RIBOTAPER_PSITES=ribotaper_psites_transcript_published_bam/P_sites_all.bed6
export HEK293_COMPARE_OUT_TAG=pub_tx
"$RSCRIPT" "$BENCH/compare_ribotaper_psite_vectors.R" . ribotaper_matched \
  2>&1 | tee compare_rt_psites_pub_tx.log
echo COMPARE_PUB_TX_DONE

# 3) Optional: matched transcript BAM if present
if [[ -f hek293_transcript_rrnafilt_primary.bam ]]; then
  if [[ ! -f ribotaper_psites_transcript_matched/P_sites_all.bed6 ]]; then
    bash "$BENCH/build_ribotaper_psites_transcript.sh" . matched_bam
  fi
  echo "=== Compare: matched transcript P-sites ==="
  export HEK293_RIBOTAPER_PSITES=ribotaper_psites_transcript_matched/P_sites_all.bed6
  export HEK293_COMPARE_OUT_TAG=matched_tx
  "$RSCRIPT" "$BENCH/compare_ribotaper_psite_vectors.R" . ribotaper_matched \
    2>&1 | tee compare_rt_psites_matched_tx.log
  echo COMPARE_MATCHED_TX_DONE
fi

# 4) Genome BED compare (sanity; uses cleaned sorted BED)
if [[ -f ribotaper_psites_published_bam/P_sites_all.bed6 ]] \
   && [[ ! -f ribotaper_psites_published_bam/P_sites_all.clean.bed6 ]]; then
  clean_genome_bed ribotaper_psites_published_bam/P_sites_all.bed6 \
    ribotaper_psites_published_bam/P_sites_all.clean.bed6
fi
if [[ -f ribotaper_psites_published_bam/P_sites_all.clean.bed6 ]]; then
  echo "=== Compare: published genome P-sites (legacy) ==="
  export HEK293_RT_PSITE_MODE=genome
  export HEK293_RIBOTAPER_PSITES=ribotaper_psites_published_bam/P_sites_all.clean.bed6
  export HEK293_COMPARE_OUT_TAG=pub_genome
  "$RSCRIPT" "$BENCH/compare_ribotaper_psite_vectors.R" . ribotaper_matched \
    2>&1 | tee compare_rt_psites_pub_genome.log || true
fi

echo "==== PUB TX ===="
tail -25 compare_rt_psites_pub_tx.log
echo ALL_COMPARE_DONE
