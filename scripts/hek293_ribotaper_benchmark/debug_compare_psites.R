#!/usr/bin/env Rscript
setwd("/home/ec2-user/work/hek293_benchmark")
rt <- read.delim("ribotaper_hek293_uorfsribotaper_matched.tsv", stringsAsFactors = FALSE)
full <- read.delim("ribostan_hek293_candidates_fullribotaper_matched.tsv",
                   stringsAsFactors = FALSE)
cat("rt cols:", paste(names(rt), collapse = ","), "\n")
cat("rt has start/start0:", "start" %in% names(rt), "start0" %in% names(rt), "\n")
cat("full p_value NA frac:", mean(is.na(full$p_value)), "\n")
if (!"start" %in% names(rt)) {
  rt$start <- rt$start0 + 1L
  rt$end <- rt$end0
}
suppressPackageStartupMessages(library(GenomicRanges))
rt_gr <- GRanges(rt$seqnames, IRanges(rt$start, rt$end), strand = rt$strand)
cand_gr <- GRanges(full$seqnames, IRanges(full$start, full$end), strand = full$strand)
hits <- findOverlaps(rt_gr, cand_gr, minoverlap = 9L)
cat("hits:", length(hits), "\n")
cmp <- read.delim("compare_ribotaper_psite_vectors_matched_genomeribotaper_matched.tsv",
                  stringsAsFactors = FALSE)
cat("cmp classes:", paste(vapply(cmp, class, ""), collapse = ","), "\n")
print(utils::head(cmp, 5))
cat("non-NA p_ribostan:", sum(!is.na(cmp$p_ribostan)), "\n")
cat("non-NA rs_orf_id:", sum(!is.na(cmp$rs_orf_id) & cmp$rs_orf_id != "NA"), "\n")
cat("C1 RT pass:", sum(cmp$fate == "C1_tested_not_significant" &
                       !is.na(cmp$p_ribotaper_psites) &
                       cmp$p_ribotaper_psites < 0.05), "\n")
cat("median n C1:", median(cmp$rt_n_psites[cmp$fate == "C1_tested_not_significant"]), "\n")
cat("median published p_sites C1:",
    median(cmp$rt_published_p_sites[cmp$fate == "C1_tested_not_significant"]), "\n")
