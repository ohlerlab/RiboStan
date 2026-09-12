#!/usr/bin/env Rscript
# Build current RiboStan uORF candidates from full GENCODE without pre-filtering.
args <- commandArgs(trailingOnly = TRUE)
workdir <- if (length(args) >= 1L) args[[1L]] else "."
setwd(workdir)

suppressPackageStartupMessages({
  library(Ribostan)
  library(GenomicRanges)
  library(rtracklayer)
})

gtf_path <- Sys.getenv("HEK293_FULL_GTF", "gencode.v19.annotation.gtf")
fasta_path <- Sys.getenv("HEK293_FASTA", "GRCh37.p13.genome.fa")
reference_path <- Sys.getenv(
  "HEK293_PARITY_UORFS",
  "ribotaper_hek293_uorfsribotaper_matched.tsv"
)
output_path <- Sys.getenv(
  "HEK293_FULL_CANDIDATES_OUTPUT",
  "ribostan_full_gencode_candidates_current.tsv"
)
start_codons <- strsplit(
  Sys.getenv("HEK293_START_CODONS", "ATG"),
  ",",
  fixed = TRUE
)[[1L]]
stopifnot(
  file.exists(gtf_path),
  file.exists(fasta_path),
  file.exists(reference_path),
  length(start_codons) > 0L,
  all(nchar(start_codons) == 3L)
)

reference <- read.delim(reference_path, stringsAsFactors = FALSE)
target_transcripts <- unique(reference$transcript_id)
gtf <- import(gtf_path, format = "gtf")
target_genes <- unique(as.character(
  gtf$gene_id[as.character(gtf$transcript_id) %in% target_transcripts]
))
keep <- as.character(gtf$transcript_id) %in% target_transcripts |
  (is.na(gtf$transcript_id) & as.character(gtf$gene_id) %in% target_genes)
gtf <- gtf[keep]
target_gtf <- tempfile(fileext = ".gtf")
on.exit(unlink(target_gtf), add = TRUE)
export(gtf, target_gtf, format = "gtf")

anno <- load_annotation(
  target_gtf,
  fasta_path,
  add_uorfs = TRUE,
  findUORFs_args = list(startCodon = start_codons)
)
candidate_ids <- names(anno$uORF)[anno$uORF]
orfs <- anno$trspacecds[candidate_ids]
stopifnot(all(width(orfs) %% 3L == 1L))

candidates <- data.frame(
  orf_id = candidate_ids,
  transcript_id = sub("_[0-9]+$", "", candidate_ids),
  tx_start = as.integer(start(orfs)),
  tx_end = as.integer(end(orfs)),
  start_codon = as.character(orfs$start_codon),
  stringsAsFactors = FALSE
)
write.table(candidates, output_path, sep = "\t", row.names = FALSE, quote = FALSE)
cat(sprintf(
  "Full-GENCODE RiboStan uORF candidates (%s): %d\n",
  paste(start_codons, collapse = ","),
  nrow(candidates)
))
message("Wrote ", output_path)
