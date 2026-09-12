#!/usr/bin/env bash
# Same as align_genome_ribotaper_style.sh, but on reads that have first been
# passed through RiboTaper's documented rRNA pre-filter (see
# filter_rrna_bowtie1.sh): GSM1887643's "Data processing" notes state "Reads
# aligning to rRNA sequences were discarded using Bowtie (v1.00), with the
# following parameters: --best" before the genome STAR alignment step. This
# was omitted from the first genome-alignment experiment; this script adds it
# back so the two pipelines are compared on equal footing.
set -euxo pipefail

WORKDIR="${1:-.}"
THREADS="${2:-30}"
cd "$WORKDIR"

FASTQ="HEK293_RIBO.norrna.fastq"
if [[ ! -f "$FASTQ" ]]; then
  echo "!!! $FASTQ not found - run filter_rrna_bowtie1.sh first" >&2
  exit 1
fi

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
  --readFilesIn "$FASTQ" \
  --outSAMtype BAM SortedByCoordinate \
  --alignEndsType EndToEnd \
  --outFilterMismatchNmax 4 \
  --outFilterMultimapNmax 8 \
  --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
  --alignSJoverhangMin 500 \
  --outSAMattributes All \
  --runThreadN "$THREADS" \
  --outFileNamePrefix hek293_genome_rrnafilt_

samtools index -@ "$THREADS" hek293_genome_rrnafilt_Aligned.sortedByCoord.out.bam
samtools view -@ "$THREADS" -b -F 0x100 \
  hek293_genome_rrnafilt_Aligned.sortedByCoord.out.bam > hek293_genome_rrnafilt_primary.bam
samtools index -@ "$THREADS" hek293_genome_rrnafilt_primary.bam

echo "=== hek293_genome_rrnafilt_Log.final.out ==="
cat hek293_genome_rrnafilt_Log.final.out
echo "=== primary-only read count ==="
samtools view -c hek293_genome_rrnafilt_primary.bam
echo "rRNA-filtered genome alignment + primary-only filter complete."
