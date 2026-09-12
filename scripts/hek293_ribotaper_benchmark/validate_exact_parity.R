#!/usr/bin/env Rscript
# Establish a zero-difference baseline on shared RiboTaper P-sites and ORFs.
args <- commandArgs(trailingOnly = TRUE)
workdir <- if (length(args) >= 1L) args[[1L]] else "."
setwd(workdir)

suppressPackageStartupMessages({
  library(Ribostan)
  library(GenomicRanges)
  library(GenomicFeatures)
  library(IRanges)
  library(multitaper)
  library(rtracklayer)
})

uorf_path <- Sys.getenv(
  "HEK293_PARITY_UORFS",
  "ribotaper_hek293_uorfsribotaper_matched.tsv"
)
psite_path <- Sys.getenv(
  "HEK293_PARITY_PSITES",
  "ribotaper_psites/P_sites_all.clean2.bed6"
)
output_path <- Sys.getenv(
  "HEK293_PARITY_OUTPUT",
  "ribotaper_ribostan_exact_parity.tsv"
)
gtf_path <- Sys.getenv("HEK293_PARITY_GTF", "gencode.v19.annotation.gtf")
stopifnot(
  file.exists(uorf_path),
  file.exists(psite_path),
  file.exists(gtf_path)
)

rt <- read.delim(uorf_path, stringsAsFactors = FALSE)
orf_ids <- unique(rt$orf_id_tr)
transcript_ids <- unique(rt$transcript_id)

message("Loading shared RiboTaper P-sites from ", psite_path, " ...")
ps <- read.delim(
  psite_path,
  header = FALSE,
  stringsAsFactors = FALSE,
  col.names = c("seqnames", "start0", "end0", "name", "score", "strand")
)
ps_gr <- GRanges(
  ps$seqnames,
  IRanges(as.integer(ps$start0) + 1L, as.integer(ps$end0)),
  strand = ps$strand
)
names(ps_gr) <- ps$name
psite_genomic_records <- data.frame(
  read_name = names(ps_gr),
  seqnames = as.character(seqnames(ps_gr)),
  position = start(ps_gr),
  strand = as.character(strand(ps_gr)),
  stringsAsFactors = FALSE
)
rm(ps)

message("Projecting the shared P-sites through RiboTaper transcript exons ...")
gtf <- import(gtf_path, format = "gtf")
exons <- gtf[
  gtf$type == "exon" &
    as.character(gtf$transcript_id) %in% transcript_ids
]
exons$transcript_id <- as.character(exons$transcript_id)
exons_by_tx <- split(exons, exons$transcript_id)
exons_by_tx <- endoapply(exons_by_tx, function(gr) {
  if (as.character(strand(gr))[1L] == "-") {
    sort(gr, decreasing = TRUE)
  } else {
    sort(gr)
  }
})
cds <- gtf[
  gtf$type == "CDS" &
    as.character(gtf$transcript_id) %in% transcript_ids
]
cds$transcript_id <- as.character(cds$transcript_id)
mapped_cds <- mapToTranscripts(cds, exons_by_tx, ignore.strand = FALSE)
same_transcript <- cds$transcript_id[mapped_cds$xHits] ==
  names(exons_by_tx)[mapped_cds$transcriptsHits]
mapped_cds <- mapped_cds[same_transcript]
cds_start_by_tx <- vapply(
  split(
    start(mapped_cds),
    names(exons_by_tx)[mapped_cds$transcriptsHits]
  ),
  min,
  integer(1L)
)
mapped <- mapToTranscripts(ps_gr, exons_by_tx, ignore.strand = FALSE)
tx_psite_records <- data.frame(
  read_name = names(ps_gr)[mapped$xHits],
  transcript_id = names(exons_by_tx)[mapped$transcriptsHits],
  position = start(mapped),
  stringsAsFactors = FALSE
)
rm(gtf, exons, cds, mapped_cds, mapped, ps_gr)
positions_by_tx <- split(
  tx_psite_records$position,
  tx_psite_records$transcript_id
)

build_window_vector <- function(transcript_id, start, end) {
  positions <- positions_by_tx[[transcript_id]]
  positions <- positions[positions >= start & positions <= end]
  as.integer(tabulate(positions - start + 1L, nbins = end - start + 1L))
}

build_vector <- function(orf_id) {
  meta <- rt[match(orf_id, rt$orf_id_tr), , drop = FALSE]
  build_window_vector(
    meta$transcript_id,
    as.integer(meta$tx_start),
    as.integer(meta$tx_end)
  )
}

vectors <- setNames(lapply(orf_ids, build_vector), orf_ids)

# Exact calculation from RiboTaper scripts/functions.R:
# take_freqs_Fvalues_all_around_3nt_spec()[c(6, 7)].
ribotaper_ftest <- function(vector, k = 24, bw = 12) {
  slepians <- dpss(
    n = if (length(vector) < 25L) 50L else length(vector),
    k = k,
    nw = bw
  )
  if (length(vector) < 25L) {
    remain <- 50L - length(vector)
    half <- as.integer(remain / 2L)
    vector <- c(rep(0, half), vector, rep(0, remain %% 2L + half))
  }
  spectrum <- spec.mtm(
    as.ts(vector),
    k = k,
    nw = bw,
    nFFT = if (length(vector) < 512L) 1024L else "default",
    centreWithSlepians = TRUE,
    Ftest = TRUE,
    maxAdaptiveIterations = 100,
    returnZeroFreq = FALSE,
    plot = FALSE,
    dpssIN = slepians
  )
  index <- which.min(abs(spectrum$freq - (1 / 3)))
  c(
    Fstat = spectrum$mtm$Ftest[index],
    spec_coef = spectrum$spec[index],
    p.value = pf(
      spectrum$mtm$Ftest[index],
      df1 = 2,
      df2 = 2 * k - 2,
      lower.tail = FALSE
    )
  )
}

run_reference <- function(vector) {
  if (sum(vector > 0L) <= 1L) {
    return(c(Fstat = NA_real_, spec_coef = NA_real_, p.value = NA_real_))
  }
  ribotaper_ftest(vector)
}
reference <- t(vapply(vectors, run_reference, numeric(3L)))

# Feed the same vectors through RiboStan's full vector-building/test path.
orf_ranges <- GRanges(
  seqnames = orf_ids,
  ranges = IRanges(1L, vapply(vectors, length, integer(1L))),
  strand = "+"
)
names(orf_ranges) <- orf_ids
anno <- list(
  trspacecds = orf_ranges,
  cdsgrl = split(orf_ranges, names(orf_ranges))
)
shared_psites <- do.call(c, lapply(orf_ids, function(orf_id) {
  vector <- vectors[[orf_id]]
  positions <- rep.int(which(vector > 0L), vector[vector > 0L])
  if (length(positions) == 0L) {
    return(GRanges())
  }
  gr <- GRanges(
    rep.int(orf_id, length(positions)),
    IRanges(positions, width = 1L),
    strand = rep.int("+", length(positions))
  )
  mcols(gr)$orf <- orf_id
  gr
}))
names(shared_psites) <- seq_along(shared_psites)

ribostan <- suppressMessages(ftest_orfs(shared_psites, anno, n_cores = 1L))
ribostan <- ribostan[match(orf_ids, ribostan$orf_id), ]

result <- data.frame(
  orf_id_tr = orf_ids,
  n_psites = vapply(vectors, sum, numeric(1L)),
  p_ribotaper = reference[, "p.value"],
  p_ribostan = ribostan$p.value,
  pct_inframe_ribotaper = vapply(vectors, function(vector) {
    if (sum(vector > 0L) <= 1L) {
      return(NA_real_)
    }
    sum(vector[seq.int(1L, length(vector), by = 3L)]) / sum(vector)
  }, numeric(1L)),
  pct_inframe_ribostan = ribostan$pct_inframe,
  stringsAsFactors = FALSE
)
result$p_difference <- result$p_ribostan - result$p_ribotaper
result$call_ribotaper <- !is.na(result$p_ribotaper) &
  result$p_ribotaper < 0.05 &
  result$pct_inframe_ribotaper > 0.5 & result$n_psites >= 6L
result$call_ribostan <- !is.na(result$p_ribostan) &
  result$p_ribostan < 0.05 &
  result$pct_inframe_ribostan > 0.5 & result$n_psites >= 6L

stopifnot(
  isTRUE(all.equal(
    result$p_ribotaper,
    result$p_ribostan,
    tolerance = 0,
    check.attributes = FALSE
  )),
  isTRUE(all.equal(
    result$pct_inframe_ribotaper,
    result$pct_inframe_ribostan,
    tolerance = 0,
    check.attributes = FALSE
  )),
  identical(result$call_ribotaper, result$call_ribostan)
)

write.table(result, output_path, sep = "\t", row.names = FALSE, quote = FALSE)
cat(sprintf(
  "EXACT PARITY: %d/%d p-values, frame fractions, and calls agree\n",
  nrow(result), nrow(result)
))
cat(sprintf("Calls on shared inputs: %d/%d\n", sum(result$call_ribostan), nrow(result)))
message("Wrote ", output_path)
