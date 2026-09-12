#!/usr/bin/env Rscript
# Stage 1: replace only published RiboTaper ORFs with RiboStan candidates.
args_all <- commandArgs(trailingOnly = FALSE)
script_path <- sub("--file=", "", args_all[grep("--file=", args_all)])
bench_dir <- if (length(script_path)) dirname(script_path) else "."
source(file.path(bench_dir, "validate_exact_parity.R"))

candidate_path <- Sys.getenv(
  "HEK293_ANNOTATION_CANDIDATES",
  "ribostan_hek293_candidates_fullribotaper_matched.tsv"
)
output_path <- Sys.getenv(
  "HEK293_ANNOTATION_DELTA_OUTPUT",
  "ribotaper_ribostan_annotation_delta.tsv"
)
stopifnot(file.exists(candidate_path))

candidates <- read.delim(candidate_path, stringsAsFactors = FALSE)
candidates <- candidates[
  !is.na(candidates$tx_start) &
    !is.na(candidates$tx_end),
  ,
  drop = FALSE
]
rt_unique <- rt[match(orf_ids, rt$orf_id_tr), , drop = FALSE]

select_candidate <- function(transcript_id, start, end) {
  rows <- candidates[candidates$transcript_id == transcript_id, , drop = FALSE]
  if (nrow(rows) == 0L) {
    return(data.frame(
      rs_orf_id = NA_character_,
      rs_tx_start = NA_integer_,
      rs_tx_end = NA_integer_,
      relation = "transcript_absent",
      stringsAsFactors = FALSE
    ))
  }

  rows$overlap <- pmax(
    0L,
    pmin(rows$tx_end, end) - pmax(rows$tx_start, start) + 1L
  )
  rows$relation <- ifelse(
    rows$tx_start == start & rows$tx_end == end,
    "exact",
    ifelse(
      rows$tx_start == start,
      "same_start",
      ifelse(rows$tx_end == end, "same_end",
             ifelse(rows$overlap > 0L, "overlap", "no_overlap"))
    )
  )
  relation_rank <- match(
    rows$relation,
    c("exact", "same_start", "same_end", "overlap", "no_overlap")
  )
  boundary_distance <- abs(rows$tx_start - start) + abs(rows$tx_end - end)
  best <- rows[order(relation_rank, boundary_distance, -rows$overlap)[1L], ]
  data.frame(
    rs_orf_id = best$orf_id,
    rs_tx_start = as.integer(best$tx_start),
    rs_tx_end = as.integer(best$tx_end),
    relation = best$relation,
    stringsAsFactors = FALSE
  )
}

matches <- do.call(rbind, Map(
  select_candidate,
  rt_unique$transcript_id,
  as.integer(rt_unique$tx_start),
  as.integer(rt_unique$tx_end)
))

evaluate_vectors <- function(input_vectors) {
  stats <- t(vapply(input_vectors, function(vector) {
    if (is.null(vector)) {
      return(c(Fstat = NA_real_, spec_coef = NA_real_, p.value = NA_real_))
    }
    run_reference(vector)
  }, numeric(3L)))
  n_psites <- vapply(input_vectors, function(vector) {
    if (is.null(vector)) 0 else sum(vector)
  }, numeric(1L))
  frame <- vapply(input_vectors, function(vector) {
    if (is.null(vector) || sum(vector > 0L) <= 1L) {
      return(NA_real_)
    }
    sum(vector[seq.int(1L, length(vector), by = 3L)]) / sum(vector)
  }, numeric(1L))
  call <- !is.na(stats[, "p.value"]) &
    stats[, "p.value"] < 0.05 &
    frame > 0.5 &
    n_psites >= 6L
  list(stats = stats, n_psites = n_psites, frame = frame, call = call)
}

# Stage 1a confirms the RiboStan window now matches RiboTaper's endpoint.
endpoint_vectors <- Map(
  build_window_vector,
  rt_unique$transcript_id,
  as.integer(rt_unique$tx_start),
  as.integer(rt_unique$tx_end)
)
endpoint <- evaluate_vectors(endpoint_vectors)

# Stage 1b then applies the generated-candidate transcript set.
candidate_vectors <- Map(function(transcript_id, start, end) {
  if (is.na(start) || is.na(end)) {
    return(NULL)
  }
  build_window_vector(transcript_id, start, end)
}, rt_unique$transcript_id, matches$rs_tx_start, matches$rs_tx_end)
candidate <- evaluate_vectors(candidate_vectors)

stage1 <- data.frame(
  orf_id_tr = orf_ids,
  transcript_id = rt_unique$transcript_id,
  rt_tx_start = as.integer(rt_unique$tx_start),
  rt_tx_end = as.integer(rt_unique$tx_end),
  baseline_n_psites = result$n_psites,
  baseline_p = result$p_ribotaper,
  baseline_pct_inframe = result$pct_inframe_ribotaper,
  baseline_call = result$call_ribotaper,
  endpoint_n_psites = endpoint$n_psites,
  endpoint_p = endpoint$stats[, "p.value"],
  endpoint_pct_inframe = endpoint$frame,
  endpoint_call = endpoint$call,
  matches,
  candidate_n_psites = candidate$n_psites,
  candidate_p = candidate$stats[, "p.value"],
  candidate_pct_inframe = candidate$frame,
  candidate_call = candidate$call,
  stringsAsFactors = FALSE
)
stage1$endpoint_outcome <- ifelse(
  stage1$baseline_call & !stage1$endpoint_call,
  "lost",
  ifelse(!stage1$baseline_call & stage1$endpoint_call, "gained", "unchanged")
)
stage1$candidate_outcome <- ifelse(
  stage1$baseline_call & !stage1$candidate_call,
  "lost",
  ifelse(!stage1$baseline_call & stage1$candidate_call, "gained", "unchanged")
)

exact <- stage1$relation == "exact"
stopifnot(
  isTRUE(all.equal(
    stage1$baseline_p[exact],
    stage1$candidate_p[exact],
    tolerance = 0,
    check.attributes = FALSE
  )),
  isTRUE(all.equal(
    stage1$baseline_pct_inframe[exact],
    stage1$candidate_pct_inframe[exact],
    tolerance = 0,
    check.attributes = FALSE
  )),
  identical(stage1$baseline_call[exact], stage1$candidate_call[exact])
)

write.table(stage1, output_path, sep = "\t", row.names = FALSE, quote = FALSE)
cat("\nEndpoint parity check:\n")
print(table(stage1$endpoint_outcome))
cat("\nSTAGE 1: replace published ORFs with selected RiboStan candidates\n")
print(table(stage1$relation))
print(table(stage1$candidate_outcome))
message("Wrote ", output_path)
