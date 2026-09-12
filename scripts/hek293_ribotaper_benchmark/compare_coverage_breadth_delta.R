#!/usr/bin/env Rscript
# Evaluate RiboTaper's 30% full-read-footprint coverage-breadth gate.
args <- commandArgs(trailingOnly = TRUE)
workdir <- if (length(args)) args[[1L]] else "."
setwd(workdir)

calls <- read.delim(
  Sys.getenv("HEK293_COVERAGE_CALLS", "ribostan_start_codon_delta.tsv"),
  stringsAsFactors = FALSE
)
qc <- read.delim(
  Sys.getenv("HEK293_COVERAGE_QC", "ribostan_multimapping_qc_delta.tsv"),
  stringsAsFactors = FALSE
)
published <- read.delim(
  Sys.getenv("HEK293_PARITY_UORFS", "ribotaper_hek293_uorfs.tsv"),
  stringsAsFactors = FALSE
)

atg <- calls[calls$scenario == "ATG" & calls$start_codon == "ATG", ]
published_key <- paste(
  published$transcript_id,
  published$tx_start,
  published$tx_end
)
atg_key <- paste(atg$transcript_id, atg$tx_start, atg$tx_end)
published_loci <- unique(atg$orf_locus_id[atg_key %in% published_key])

called <- atg[atg$locus_representative & atg$call_bh, ]
called <- merge(
  called,
  qc[c("orf_id", "fraction_all", "fraction_unique", "only_multi_ratio")],
  by = "orf_id"
)
called$width <- called$tx_end - called$tx_start + 1L
called$published <- called$orf_locus_id %in% published_loci
called$passes_30pct <- called$fraction_all > 0.3

thresholds <- seq(0, 1, by = 0.05)
threshold_summary <- do.call(rbind, lapply(thresholds, function(threshold) {
  retained <- called$fraction_all > threshold
  data.frame(
    threshold = threshold,
    retained = sum(retained),
    removed = sum(!retained),
    published_retained = sum(called$published & retained),
    published_removed = sum(called$published & !retained),
    nonpublished_retained = sum(!called$published & retained),
    nonpublished_removed = sum(!called$published & !retained)
  )
}))

output_path <- Sys.getenv(
  "HEK293_COVERAGE_OUTPUT",
  "ribostan_coverage_breadth_delta.tsv"
)
summary_path <- Sys.getenv(
  "HEK293_COVERAGE_SUMMARY",
  "ribostan_coverage_breadth_thresholds.tsv"
)
write.table(called, output_path, sep = "\t", row.names = FALSE, quote = FALSE)
write.table(
  threshold_summary,
  summary_path,
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

print(threshold_summary)
print(with(called, table(published, passes_30pct)))
print(with(
  called,
  aggregate(
    cbind(width, n_psites, p.value, q.value, pct_inframe),
    list(passes_30pct = passes_30pct),
    median
  )
))
comparison <- with(called, table(published, passes_30pct))
if (all(dim(comparison) >= 2L)) {
  print(fisher.test(comparison))
}
message("Wrote ", output_path, " and ", summary_path)
