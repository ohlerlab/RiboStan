#!/usr/bin/env bash
# Align SRR2433794 to the genome the way RiboTaper's own HEK293 analysis did
# (GSM1887643 "Data processing" notes): STAR against hg19/GRCh37 with
# Gencode v19 splice junctions, then keep only the primary alignment per read
# (samtools -F 0x100), matching "Alignments flagged as secondary alignment
# were filtered out ... ensuring 1 genomic position per aligned read."
#
# This produces a genome-space BAM; a separate R step
# (convert_genome_to_txspace.R) projects it into the transcript-space
# coordinates Ribostan's periodicity pipeline expects.
set -euxo pipefail

WORKDIR="${1:-.}"
THREADS="${2:-30}"
cd "$WORKDIR"

if [[ ! -f star_genome_index/SAindex ]]; then
  mkdir -p star_genome_index
  STAR --runMode genomeGenerate \
    --genomeDir star_genome_index \
    --genomeFastaFiles GRCh37.p13.genome.fa \
    --sjdbGTFfile gencode.v19.annotation.gtf \
    --sjdbOverhang 29 \
    --runThreadN "$THREADS"
fi

STAR --genomeDir star_genome_index \
  --readFilesIn HEK293_RIBO.fastq.gz \
  --readFilesCommand zcat \
  --outSAMtype BAM SortedByCoordinate \
  --alignEndsType EndToEnd \
  --outFilterMismatchNmax 4 \
  --outFilterMultimapNmax 8 \
  --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
  --alignSJoverhangMin 500 \
  --outSAMattributes All \
  --runThreadN "$THREADS" \
  --outFileNamePrefix hek293_genome_

samtools index -@ "$THREADS" hek293_genome_Aligned.sortedByCoord.out.bam
samtools view -@ "$THREADS" -b -F 0x100 \
  hek293_genome_Aligned.sortedByCoord.out.bam > hek293_genome_primary.bam
samtools index -@ "$THREADS" hek293_genome_primary.bam

echo "=== hek293_genome_Log.final.out ==="
cat hek293_genome_Log.final.out
echo "=== primary-only read count ==="
samtools view -c hek293_genome_primary.bam
echo "Genome alignment + primary-only filter complete."
