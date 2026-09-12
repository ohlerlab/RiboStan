#!/usr/bin/env bash
# STAR transcriptome alignment for HEK293 Ribo-seq (Ribostan expects transcript-space BAM).
set -euxo pipefail

WORKDIR="${1:-.}"
N_CORES="${2:-8}"
cd "$WORKDIR"

STAR_BIN="${STAR_BIN:-STAR}"
SAMTOOLS="${SAMTOOLS:-samtools}"

STAR_RAM="${STAR_RAM:-80000000000}"

if [[ ! -d star_tx_index ]] || [[ ! -f star_tx_index/genomeParameters.txt ]]; then
  rm -rf star_tx_index
  mkdir -p star_tx_index
  "$STAR_BIN" --runMode genomeGenerate \
    --genomeDir star_tx_index \
    --genomeFastaFiles gencode.v19.transcripts.fa \
    --runThreadN "$N_CORES" \
    --genomeSAindexNbases 12 \
    --limitGenomeGenerateRAM "$STAR_RAM"
fi

"$STAR_BIN" --runMode alignReads \
  --genomeDir star_tx_index \
  --readFilesIn HEK293_RIBO.fastq.gz \
  --readFilesCommand zcat \
  --outSAMtype BAM SortedByCoordinate \
  --outFilterMultimapNmax 20 \
  --outFilterMismatchNmax 4 \
  --alignEndsType EndToEnd \
  --runThreadN "$N_CORES" \
  --outFileNamePrefix hek293_tx_ \
  --limitBAMsortRAM 8000000000

"$SAMTOOLS" index hek293_tx_Aligned.sortedByCoord.out.bam
ln -sf hek293_tx_Aligned.sortedByCoord.out.bam HEK293_ribo_transcriptome.bam
ln -sf hek293_tx_Aligned.sortedByCoord.out.bam.bai HEK293_ribo_transcriptome.bam.bai

echo "Transcriptome BAM: $(pwd)/HEK293_ribo_transcriptome.bam"
