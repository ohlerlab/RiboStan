#!/usr/bin/env Rscript
# Parse RiboTaper HEK293 uORF calls from GEO supplementary BED (GSE73136).
#
# Calviello et al. 2016 (Nature Methods): >600 genes with translated uORFs;
# GSM1887643_HEK293_filtered_ORFs.bed.gz contains periodic ORFs including 908
# uORF rows (701 unique transcripts).
load_ribotaper_hek293_uorfs <- function(bed_gz) {
  if (!file.exists(bed_gz)) {
    stop("Missing RiboTaper BED: ", bed_gz)
  }
  con <- gzfile(bed_gz, open = "rt")
  on.exit(close(con), add = TRUE)
  raw <- read.delim(con, header = FALSE, stringsAsFactors = FALSE,
                    comment.char = "", quote = "")

  # BED6+: chr, start (0-based), end (exclusive), name, score, strand
  names(raw) <- c("seqnames", "start0", "end0", "name", "p_sites", "strand")[seq_len(ncol(raw))]

  meta <- strsplit(raw$name, ";", fixed = TRUE)
  raw$orf_id_tr <- vapply(meta, `[`, character(1), 1)
  raw$category <- vapply(meta, function(x) if (length(x) >= 2) x[[2]] else NA_character_, character(1))
  raw$gene_biotype <- vapply(meta, function(x) if (length(x) >= 3) x[[3]] else NA_character_, character(1))

  parts <- strsplit(raw$orf_id_tr, "_", fixed = TRUE)
  raw$transcript_id <- vapply(parts, `[`, character(1), 1)
  raw$tx_start <- vapply(parts, function(p) {
    if (length(p) < 3) return(NA_integer_)
    as.integer(p[[length(p) - 1]])
  }, integer(1))
  raw$tx_end <- vapply(parts, function(p) {
    if (length(p) < 3) return(NA_integer_)
    as.integer(p[[length(p)]])
  }, integer(1))

  # Convert BED to 1-based closed GRanges coordinates.
  raw$start <- as.integer(raw$start0) + 1L
  raw$end <- as.integer(raw$end0)
  raw$width_nt <- raw$end - raw$start + 1L
  raw$width_codons <- (raw$width_nt - 1L) / 3L

  uorfs <- raw[raw$category == "uORF", , drop = FALSE]
  rownames(uorfs) <- NULL
  uorfs
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  workdir <- if (length(args) >= 1) args[[1]] else "."
  bed <- file.path(workdir, "GSM1887643_HEK293_filtered_ORFs.bed.gz")
  uorfs <- load_ribotaper_hek293_uorfs(bed)
  out_tsv <- file.path(workdir, "ribotaper_hek293_uorfs.tsv")
  write.table(uorfs, out_tsv, sep = "\t", row.names = FALSE, quote = FALSE)
  saveRDS(uorfs, file.path(workdir, "ribotaper_hek293_uorfs.rds"))
  message("RiboTaper HEK293 uORFs: ", nrow(uorfs), " (unique tx: ",
          length(unique(uorfs$transcript_id)), ")")
}
