#!/usr/bin/env Rscript
# Isolate Ribostan's CDS-calibrated phase reranking on fixed P-site records.
args <- commandArgs(trailingOnly = TRUE)
workdir <- if (length(args)) args[[1L]] else "."
setwd(workdir)

suppressPackageStartupMessages({
  library(GenomicFeatures)
  library(GenomicRanges)
  library(Ribostan)
  library(rtracklayer)
})

orf_path <- Sys.getenv(
  "HEK293_PARITY_UORFS",
  "ribotaper_hek293_uorfsribotaper_matched.tsv"
)
record_path <- Sys.getenv(
  "HEK293_PHASE_RECORDS",
  "ribotaper_ribostan_psite_record_delta_stage4_lengths.tsv"
)
gtf_path <- Sys.getenv("HEK293_PARITY_GTF", "gencode.v19.annotation.gtf")
output_path <- Sys.getenv(
  "HEK293_PHASE_OUTPUT",
  "ribostan_phase_shift_delta.tsv"
)
calibration_path <- Sys.getenv(
  "HEK293_PHASE_CALIBRATION_OUTPUT",
  "ribostan_phase_shift_calibration.tsv"
)

orfs <- read.delim(orf_path, stringsAsFactors = FALSE)
orfs <- orfs[!duplicated(orfs$orf_id_tr), , drop = FALSE]
records <- read.delim(record_path, stringsAsFactors = FALSE)
records <- records[
  !is.na(records$position_ribostan) & records$read_length %in% 26:30,
  ,
  drop = FALSE
]

gtf <- import(gtf_path, format = "gtf")
transcript_ids <- unique(orfs$transcript_id)
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
mapped_cds_tx <- names(exons_by_tx)[mapped_cds$transcriptsHits]
cds_starts <- aggregate(
  start(mapped_cds),
  by = list(transcript_id = mapped_cds_tx),
  FUN = min
)
names(cds_starts)[[2L]] <- "cds_start"
cds_ends <- aggregate(
  end(mapped_cds),
  by = list(transcript_id = mapped_cds_tx),
  FUN = max
)
names(cds_ends)[[2L]] <- "cds_end"
cds_bounds <- merge(cds_starts, cds_ends, by = "transcript_id")

records$cds_start <- NULL
records <- merge(records, cds_bounds, by = "transcript_id", all.x = TRUE)
records$read_end <- records$read_start + records$read_length - 1L
calibration <- records[
  !is.na(records$cds_start) &
    records$read_start <= records$cds_end &
    records$read_end >= records$cds_start,
  ,
  drop = FALSE
]
calibration$phase <- (
  calibration$read_start - calibration$cds_start
) %% 3L
shift_df <- aggregate(
  list(n = rep.int(1L, nrow(calibration))),
  by = list(
    read_length = calibration$read_length,
    phase = calibration$phase
  ),
  FUN = sum
)
shift_df <- do.call(rbind, lapply(
  split(shift_df, shift_df$read_length),
  function(group) {
    group$shift <- rank(-group$n) - 1L
    group
  }
))
stopifnot(all(shift_df$shift %in% 0:2))
write.table(
  shift_df,
  calibration_path,
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

records_by_tx <- split(records, records$transcript_id)
evaluate_orf <- function(orf, apply_shift) {
  selected <- records_by_tx[[orf$transcript_id]]
  selected <- selected[
    selected$read_start <= orf$tx_end &
      selected$read_end >= orf$tx_start,
    ,
    drop = FALSE
  ]
  positions <- selected$position_ribostan
  if (apply_shift && nrow(selected)) {
    phase <- (selected$read_start - orf$tx_start) %% 3L
    shift_index <- match(
      paste(selected$read_length, phase),
      paste(shift_df$read_length, shift_df$phase)
    )
    stopifnot(!anyNA(shift_index))
    positions <- positions - phase + shift_df$shift[shift_index]
  }
  positions <- positions[positions >= orf$tx_start & positions <= orf$tx_end]
  vector <- as.integer(tabulate(
    positions - orf$tx_start + 1L,
    nbins = orf$tx_end - orf$tx_start + 1L
  ))
  n_psites <- sum(vector)
  if (sum(vector > 0L) <= 1L) {
    return(c(n_psites = n_psites, p = NA_real_, frame = NA_real_, call = 0))
  }
  p <- unname(Ribostan:::ftestvect(vector)["p.value"])
  frame <- sum(vector[seq.int(1L, length(vector), by = 3L)]) / n_psites
  call <- p < 0.05 && frame > 0.5 && n_psites >= 6L
  c(n_psites = n_psites, p = p, frame = frame, call = call)
}

evaluate_mode <- function(apply_shift) {
  result <- t(vapply(
    seq_len(nrow(orfs)),
    function(index) evaluate_orf(orfs[index, , drop = FALSE], apply_shift),
    numeric(4L)
  ))
  rownames(result) <- orfs$orf_id_tr
  result
}

unshifted <- evaluate_mode(FALSE)
shifted <- evaluate_mode(TRUE)
unshifted_call <- as.logical(unshifted[, "call"])
shifted_call <- as.logical(shifted[, "call"])

prior_path <- "ribostan_read_length_subset_delta.tsv"
if (file.exists(prior_path)) {
  prior <- read.delim(prior_path, stringsAsFactors = FALSE)
  prior <- prior[prior$subset == "add_30", ]
  prior <- prior[match(rownames(unshifted), prior$orf_id_tr), ]
  stopifnot(identical(unshifted_call, prior$call))
}

output <- data.frame(
  orf_id_tr = rownames(unshifted),
  n_unshifted = unshifted[, "n_psites"],
  n_shifted = shifted[, "n_psites"],
  p_unshifted = unshifted[, "p"],
  p_shifted = shifted[, "p"],
  frame_unshifted = unshifted[, "frame"],
  frame_shifted = shifted[, "frame"],
  call_unshifted = unshifted_call,
  call_shifted = shifted_call,
  stringsAsFactors = FALSE
)
output$outcome <- ifelse(
  output$call_unshifted & !output$call_shifted,
  "lost",
  ifelse(
    !output$call_unshifted & output$call_shifted,
    "gained",
    "unchanged"
  )
)
write.table(
  output,
  output_path,
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

print(shift_df[order(shift_df$read_length, shift_df$phase), ], row.names = FALSE)
print(table(output$outcome))
message("Wrote ", output_path)
