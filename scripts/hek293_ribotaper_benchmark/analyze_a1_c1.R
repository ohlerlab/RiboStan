#!/usr/bin/env Rscript
# Summarize remaining A1 (transcript absent) and C1 (borderline p) gaps.
args <- commandArgs(trailingOnly = TRUE)
workdir <- if (length(args) >= 1) args[[1]] else "."
suffix <- if (length(args) >= 2) args[[2]] else Sys.getenv("HEK293_OUT_SUFFIX", "ribotaper_matched")
setwd(workdir)

p <- function(base) sub("\\.", paste0(suffix, "."), paste0(base, ".tsv"))

fate <- read.delim(sub("\\.", paste0(suffix, "."), "diag_ribotaper_fate.tsv"),
                   stringsAsFactors = FALSE)
cand <- read.delim(p("ribostan_hek293_candidates_full"), stringsAsFactors = FALSE)
at <- read.delim(p("ribostan_hek293_all_tested"), stringsAsFactors = FALSE)
at_p <- at[, c("orf_id", "p_value", "n_psites", "pct_inframe")]
full <- merge(cand, at_p, by = "orf_id", all.x = TRUE)

rt_gr <- fate[fate$class == "C1_tested_not_significant", , drop = FALSE]
c1_orfs <- full[full$transcript_id %in% rt_gr$transcript_id & !is.na(full$p_value), ]
c1_orfs <- c1_orfs[order(c1_orfs$transcript_id, c1_orfs$p_value), ]
c1_best <- c1_orfs[!duplicated(c1_orfs$transcript_id), ]

cat("=== C1_tested_not_significant (n=", nrow(rt_gr), " uORFs, ",
    length(unique(rt_gr$transcript_id)), " transcripts) ===\n", sep = "")
cat(sprintf("Best overlapping candidate per C1 uORF:\n"))
cat(sprintf("  median p = %.4g (RECOVERED median = %.4g)\n",
            median(c1_best$p_value, na.rm = TRUE),
            median(full$p_value[full$orf_id %in% read.delim(p("ribostan_hek293_periodic"))$orf_id],
                   na.rm = TRUE)))
cat(sprintf("  p in [0.05, 0.10): %d\n",
            sum(c1_best$p_value >= 0.05 & c1_best$p_value < 0.10, na.rm = TRUE)))
cat(sprintf("  p in [0.10, 0.20): %d\n",
            sum(c1_best$p_value >= 0.10 & c1_best$p_value < 0.20, na.rm = TRUE)))
cat(sprintf("  p >= 0.20: %d\n",
            sum(c1_best$p_value >= 0.20, na.rm = TRUE)))
cat(sprintf("  median n_psites = %.0f, median pct_inframe = %.2f\n",
            median(c1_best$n_psites, na.rm = TRUE),
            median(c1_best$pct_inframe, na.rm = TRUE)))

a1_tx <- unique(fate$transcript_id[fate$class == "A1_transcript_absent"])
wl <- readLines("ccds_appris_transcript_whitelist.txt")
cat("\n=== A1_transcript_absent (n=", sum(fate$class == "A1_transcript_absent"),
    " uORFs, ", length(a1_tx), " transcripts) ===\n", sep = "")
cat(sprintf("In filtered whitelist: %d / %d\n", sum(a1_tx %in% wl), length(a1_tx)))

gtf_lines <- readLines("gencode.v19.annotation.gtf", warn = FALSE)
tx_line <- function(tx) {
  gtf_lines[grepl(paste0('transcript_id "', tx, '"'), gtf_lines) &
              grepl("\ttranscript\t", gtf_lines)][1]
}

rows <- lapply(a1_tx, function(tx) {
  line <- tx_line(tx)
  gid <- if (is.na(line)) NA_character_ else {
    sub('.*gene_id "([^"]+)".*', "\\1", line)
  }
  prin_tx <- if (is.na(gid)) character(0) else {
    gl <- gtf_lines[grepl(paste0('gene_id "', gid, '"'), gtf_lines)]
    prin <- gl[grepl("\ttranscript\t", gl) & grepl("appris_prin", gl)]
    unique(sub('.*transcript_id "([^"]+)".*', "\\1", prin))
  }
  kept_prin <- intersect(prin_tx, wl)
  data.frame(
    transcript_id = tx,
    in_whitelist = tx %in% wl,
    gene_id = gid,
    appris_principal = paste(prin_tx, collapse = ";"),
    kept_principal = paste(kept_prin, collapse = ";"),
    has_ccds = !is.na(line) && grepl("ccdsid", line, ignore.case = TRUE),
    stringsAsFactors = FALSE
  )
})
a1_df <- do.call(rbind, rows)
print(a1_df)
write.table(a1_df, sub("\\.", paste0(suffix, "."), "diag_a1_transcript_detail.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
