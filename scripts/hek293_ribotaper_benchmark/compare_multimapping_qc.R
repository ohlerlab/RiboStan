#!/usr/bin/env Rscript
# Apply RiboTaper's genomic read-footprint multimapping QC to fixed calls.
args <- commandArgs(trailingOnly = TRUE)
workdir <- if (length(args)) args[[1L]] else "."
setwd(workdir)

suppressPackageStartupMessages({
  library(Ribostan)
  library(GenomicRanges)
  library(rtracklayer)
})

reference <- read.delim(
  Sys.getenv("HEK293_PARITY_UORFS", "ribotaper_hek293_uorfs.tsv"),
  stringsAsFactors = FALSE
)
gtf_path <- Sys.getenv("HEK293_FULL_GTF", "gencode.v19.annotation.gtf")
fasta_path <- Sys.getenv("HEK293_FASTA", "GRCh37.p13.genome.fa")
bam_path <- Sys.getenv(
  "HEK293_GENOME_BAM",
  "hek293_genome_rrnafilt_primary.bam"
)
call_path <- Sys.getenv(
  "HEK293_MULTIMAP_CALLS",
  "ribostan_start_codon_delta.tsv"
)
output_path <- Sys.getenv(
  "HEK293_MULTIMAP_OUTPUT",
  "ribostan_multimapping_qc_delta.tsv"
)
bedtools <- Sys.getenv(
  "BEDTOOLS",
  "/opt/micromamba/envs/rnaseq/bin/bedtools"
)
samtools <- Sys.getenv(
  "SAMTOOLS",
  "/opt/micromamba/envs/rnaseq/bin/samtools"
)
n_cores <- Sys.getenv("HEK293_N_CORES", "8")

target_transcripts <- unique(reference$transcript_id)
gtf <- import(gtf_path, format = "gtf")
target_genes <- unique(as.character(
  gtf$gene_id[as.character(gtf$transcript_id) %in% target_transcripts]
))
keep <- as.character(gtf$transcript_id) %in% target_transcripts |
  (is.na(gtf$transcript_id) & as.character(gtf$gene_id) %in% target_genes)
target_gtf <- tempfile(fileext = ".gtf")
on.exit(unlink(target_gtf), add = TRUE)
export(gtf[keep], target_gtf, format = "gtf")
rm(gtf)

anno <- load_annotation(target_gtf, fasta_path, add_uorfs = TRUE)
uorf_ids <- names(anno$uORF)[anno$uORF]
genomic_orfs <- anno$cdsgrl[uorf_ids]
bed_path <- tempfile(fileext = ".bed")
on.exit(unlink(bed_path), add = TRUE)
bed_rows <- lapply(seq_along(genomic_orfs), function(index) {
  blocks <- sort(genomic_orfs[[index]])
  data.frame(
    chrom = as.character(seqnames(blocks)),
    chrom_start = start(blocks) - 1L,
    chrom_end = end(blocks),
    name = names(genomic_orfs)[index],
    score = 0L,
    strand = as.character(strand(blocks))
  )
})
write.table(
  do.call(rbind, bed_rows),
  bed_path,
  sep = "\t",
  row.names = FALSE,
  col.names = FALSE,
  quote = FALSE
)

unique_bam <- "hek293_genome_rrnafilt_unique_mapq50.bam"
if (!file.exists(unique_bam)) {
  status <- system2(
    samtools,
    c("view", "-@", n_cores, "-b", "-q", "50", "-o", unique_bam, bam_path)
  )
  stopifnot(status == 0L)
}

run_coverage <- function(bam, output) {
  status <- system2(
    bedtools,
    c("coverage", "-a", bed_path, "-b", bam, "-s", "-split"),
    stdout = output
  )
  stopifnot(status == 0L)
  coverage <- read.delim(output, header = FALSE, stringsAsFactors = FALSE)
  block_coverage <- data.frame(
    orf_id = coverage[[4L]],
    covered = coverage[[ncol(coverage) - 2L]],
    width = coverage[[ncol(coverage) - 1L]],
    stringsAsFactors = FALSE
  )
  totals <- aggregate(cbind(covered, width) ~ orf_id, block_coverage, sum)
  data.frame(
    orf_id = totals$orf_id,
    fraction = totals$covered / totals$width,
    stringsAsFactors = FALSE
  )
}

all_coverage <- run_coverage(bam_path, tempfile(fileext = ".tsv"))
unique_coverage <- run_coverage(unique_bam, tempfile(fileext = ".tsv"))
qc <- merge(
  all_coverage,
  unique_coverage,
  by = "orf_id",
  suffixes = c("_all", "_unique")
)
qc$fraction_only_multi <- qc$fraction_all - qc$fraction_unique
qc$only_multi_ratio <- ifelse(
  qc$fraction_all > 0,
  qc$fraction_only_multi / qc$fraction_all,
  NA_real_
)
qc$low_coverage <- qc$fraction_all <= 0.3
qc$multimap_dependent <- !is.na(qc$only_multi_ratio) &
  qc$only_multi_ratio >= 0.3
qc$flagged <- qc$low_coverage | qc$multimap_dependent

calls <- read.delim(call_path, stringsAsFactors = FALSE)
calls <- calls[
  calls$scenario == "ATG" &
    calls$start_codon == "ATG" &
    calls$locus_representative,
  ,
  drop = FALSE
]
qc <- merge(
  qc,
  calls[c("orf_id", "orf_locus_id", "call_bh")],
  by = "orf_id",
  all.x = TRUE
)
write.table(qc, output_path, sep = "\t", row.names = FALSE, quote = FALSE)

called <- qc[!is.na(qc$call_bh) & qc$call_bh, , drop = FALSE]
print(table(
  low_coverage = called$low_coverage,
  multimap_dependent = called$multimap_dependent
))
cat(sprintf(
  "BH calls retained after QC: %d/%d\n",
  sum(!called$flagged),
  nrow(called)
))
message("Wrote ", output_path)
