#!/usr/bin/env Rscript
# Isolate the effects of adding 25 nt and 30 nt reads to the 26-29 nt baseline.
args <- commandArgs(trailingOnly = TRUE)
workdir <- if (length(args)) args[[1L]] else "."
setwd(workdir)

suppressPackageStartupMessages(library(Ribostan))

orf_path <- Sys.getenv(
  "HEK293_PARITY_UORFS",
  "ribotaper_hek293_uorfsribotaper_matched.tsv"
)
record_path <- Sys.getenv(
  "HEK293_LENGTH_RECORDS",
  "ribotaper_ribostan_psite_record_delta_stage4_lengths.tsv"
)
output_path <- Sys.getenv(
  "HEK293_LENGTH_SUBSET_OUTPUT",
  "ribostan_read_length_subset_delta.tsv"
)

orfs <- read.delim(orf_path, stringsAsFactors = FALSE)
orfs <- orfs[!duplicated(orfs$orf_id_tr), , drop = FALSE]
records <- read.delim(record_path, stringsAsFactors = FALSE)
records <- records[!is.na(records$position_ribostan), , drop = FALSE]

evaluate_subset <- function(read_lengths) {
  selected <- records[records$read_length %in% read_lengths, , drop = FALSE]
  positions_by_tx <- split(selected$position_ribostan, selected$transcript_id)
  result <- t(vapply(seq_len(nrow(orfs)), function(index) {
    orf <- orfs[index, , drop = FALSE]
    positions <- positions_by_tx[[orf$transcript_id]]
    positions <- positions[
      positions >= orf$tx_start & positions <= orf$tx_end
    ]
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
  }, numeric(4L)))
  rownames(result) <- orfs$orf_id_tr
  result
}

subsets <- list(
  baseline_26_29 = 26:29,
  add_25 = 25:29,
  add_30 = 26:30,
  add_25_30 = 25:30
)
results <- lapply(subsets, evaluate_subset)
baseline_call <- as.logical(results$baseline_26_29[, "call"])

output <- do.call(rbind, lapply(names(results), function(subset_name) {
  result <- results[[subset_name]]
  call <- as.logical(result[, "call"])
  data.frame(
    orf_id_tr = rownames(result),
    subset = subset_name,
    n_psites = result[, "n_psites"],
    p = result[, "p"],
    frame = result[, "frame"],
    call = call,
    outcome_vs_baseline = ifelse(
      baseline_call & !call,
      "lost",
      ifelse(!baseline_call & call, "gained", "unchanged")
    ),
    stringsAsFactors = FALSE
  )
}))
write.table(
  output,
  output_path,
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

print(with(
  output[output$subset != "baseline_26_29", ],
  table(subset, outcome_vs_baseline)
))
message("Wrote ", output_path)
