#!/usr/bin/env Rscript
# Stage 2: compare P-site construction with reads, lengths, and offsets fixed.
args_all <- commandArgs(trailingOnly = FALSE)
script_path <- sub("--file=", "", args_all[grep("--file=", args_all)])
bench_dir <- if (length(script_path)) dirname(script_path) else "."
source(file.path(bench_dir, "validate_exact_parity.R"))

bam_path <- Sys.getenv(
  "HEK293_GENOME_BAM",
  "hek293_genome_rrnafilt_primary.bam"
)
record_output <- Sys.getenv(
  "HEK293_PSITE_RECORD_OUTPUT",
  "ribotaper_ribostan_psite_record_delta.tsv"
)
orf_output <- Sys.getenv(
  "HEK293_PSITE_ORF_OUTPUT",
  "ribotaper_ribostan_psite_orf_delta.tsv"
)
stopifnot(file.exists(bam_path))

suppressPackageStartupMessages({
  library(GenomicAlignments)
  library(Rsamtools)
})
gal <- readGAlignments(
  bam_path,
  use.names = TRUE,
  param = ScanBamParam(flag = scanBamFlag(isUnmappedQuery = FALSE))
)
match_baseline_reads <- identical(
  toupper(Sys.getenv("HEK293_PSITE_MATCH_BASELINE_READS", "TRUE")),
  "TRUE"
)
if (match_baseline_reads) {
  selected_read_names <- unique(psite_genomic_records$read_name)
  gal <- gal[names(gal) %in% selected_read_names]
}
cigar_ops <- explodeCigarOps(cigar(gal))
cigar_lengths <- explodeCigarOpLengths(cigar(gal))
ribotaper_compatible <- identical(
  toupper(Sys.getenv("HEK293_PSITE_RIBOTAPER_COMPATIBLE", "TRUE")),
  "TRUE"
)

get_bed_blocks <- function(ops, lengths) {
  block_starts <- 0L
  block_sizes <- 0L
  reference_cursor <- 0L
  block_index <- 1L
  for (index in seq_along(ops)) {
    op <- ops[[index]]
    op_length <- lengths[[index]]
    if (op == "N") {
      reference_cursor <- reference_cursor + op_length
      block_index <- block_index + 1L
      block_starts[[block_index]] <- reference_cursor
      block_sizes[[block_index]] <- 0L
    } else if (op %in% c("M", "D", "=", "X")) {
      block_sizes[[block_index]] <- block_sizes[[block_index]] + op_length
      reference_cursor <- reference_cursor + op_length
    }
  }
  list(starts = block_starts, sizes = block_sizes)
}
bed_blocks <- mapply(
  get_bed_blocks,
  cigar_ops,
  cigar_lengths,
  SIMPLIFY = FALSE
)
selection_width <- vapply(
  bed_blocks,
  function(blocks) {
    sizes <- blocks$sizes[seq_len(min(3L, length(blocks$sizes)))]
    sum(sizes)
  },
  integer(1L)
)
selected_lengths <- as.integer(strsplit(
  Sys.getenv("HEK293_PSITE_LENGTHS", "26,27,28,29"),
  ",",
  fixed = TRUE
)[[1L]])
stopifnot(!anyNA(selected_lengths))
keep_length <- selection_width %in% selected_lengths
gal <- gal[keep_length]
cigar_ops <- cigar_ops[keep_length]
cigar_lengths <- cigar_lengths[keep_length]
bed_blocks <- bed_blocks[keep_length]
selection_width <- selection_width[keep_length]

offset_path <- Sys.getenv("HEK293_PSITE_OFFSETS", "")
if (nzchar(offset_path)) {
  offsets <- read.table(offset_path, header = TRUE, stringsAsFactors = FALSE)
  offsets <- unique(offsets[c("readlen", "p_offset")])
  stopifnot(!anyDuplicated(offsets$readlen))
  offset_by_record <- offsets$p_offset[match(selection_width, offsets$readlen)]
  stopifnot(!anyNA(offset_by_record))
} else {
  offset_by_record <- rep.int(12L, length(gal))
}

map_aligned_offset <- function(reference_start, strand, ops, lengths, offset) {
  reference_cursor <- reference_start
  indices <- if (strand == "+") seq_along(ops) else rev(seq_along(ops))
  if (strand == "-") {
    reference_cursor <- reference_start +
      sum(lengths[ops %in% c("M", "D", "N", "=", "X")]) - 1L
  }
  for (index in indices) {
    op <- ops[[index]]
    op_length <- lengths[[index]]
    aligned_reference <- op %in% c("M", "D", "=", "X")
    if (aligned_reference && offset < op_length) {
      direction <- if (strand == "+") 1L else -1L
      return(reference_cursor + direction * offset)
    }
    if (aligned_reference) {
      offset <- offset - op_length
    }
    if (op %in% c("M", "D", "N", "=", "X")) {
      direction <- if (strand == "+") 1L else -1L
      reference_cursor <- reference_cursor + direction * op_length
    }
  }
  NA_integer_
}

map_ribotaper_offset <- function(reference_start, strand, blocks) {
  sizes <- c(blocks$sizes, rep.int(0L, 4L))[seq_len(4L)]
  starts <- c(blocks$starts, rep.int(0L, 4L))[seq_len(4L)]
  total <- sum(blocks$sizes)
  position <- reference_start - 1L
  if (strand == "+") {
    if (sizes[[1L]] < 12L && sum(sizes[1:2]) > 12L) {
      position <- position + starts[[2L]] + 12L - sizes[[1L]]
    }
    if (sizes[[1L]] >= 12L) {
      position <- position + 12L
    }
    if (sizes[[1L]] < 12L &&
        sum(sizes[1:2]) < 12L &&
        sum(sizes[1:3]) > 12L) {
      position <- position + starts[[3L]] + 12L - sum(sizes[1:2])
    }
    if (sum(sizes[1:3]) < 12L) {
      position <- position + starts[[4L]] + 12L - sum(sizes[1:3])
    }
  } else {
    target <- total - 12L
    if (sizes[[1L]] >= target) {
      position <- position + target - 1L
    }
    if (sizes[[1L]] < target && sum(sizes[1:2]) >= target) {
      position <- position + starts[[2L]] + target - 1L - sizes[[1L]]
    }
    if (sum(sizes[1:2]) < target && sum(sizes[1:3]) >= target) {
      position <- position + starts[[3L]] + target - 1L - sum(sizes[1:2])
    }
  }
  position + 1L
}
if (ribotaper_compatible) {
  reference_position <- mapply(
    map_ribotaper_offset,
    start(gal),
    as.character(strand(gal)),
    bed_blocks,
    USE.NAMES = FALSE
  )
} else {
  reference_position <- mapply(
    map_aligned_offset,
    start(gal),
    as.character(strand(gal)),
    cigar_ops,
    cigar_lengths,
    offset_by_record,
    USE.NAMES = FALSE
  )
}
keep <- !is.na(reference_position)
rs_psites <- GRanges(
  seqnames(gal)[keep],
  IRanges(reference_position[keep], width = 1L),
  strand = strand(gal)[keep]
)
names(rs_psites) <- names(gal)[keep]
seqlevels(rs_psites) <- sub("\\|.*", "", seqlevels(rs_psites))
psite_read_length <- selection_width[keep]
psite_offset <- offset_by_record[keep]
rs_mapped <- mapToTranscripts(
  rs_psites,
  exons_by_tx,
  ignore.strand = FALSE
)
rs_records <- data.frame(
  read_name = names(rs_psites)[rs_mapped$xHits],
  transcript_id = names(exons_by_tx)[rs_mapped$transcriptsHits],
  position = as.integer(start(rs_mapped)),
  read_length = psite_read_length[rs_mapped$xHits],
  psite_offset = psite_offset[rs_mapped$xHits],
  stringsAsFactors = FALSE
)

rs_records$read_start <- rs_records$position - rs_records$psite_offset
rs_records$cds_start <- unname(
  cds_start_by_tx[rs_records$transcript_id]
)
rs_records$start_distance <- rs_records$read_start - rs_records$cds_start
calibration_records <- rs_records[
  !is.na(rs_records$start_distance) &
    rs_records$start_distance >= -35L &
    rs_records$start_distance <= -5L,
  ,
  drop = FALSE
]
offset_histogram <- aggregate(
  list(n = rep.int(1L, nrow(calibration_records))),
  by = list(
    readlen = calibration_records$read_length,
    distance = calibration_records$start_distance
  ),
  FUN = sum
)
offset_histogram <- offset_histogram[
  order(
    offset_histogram$readlen,
    -offset_histogram$n,
    offset_histogram$distance
  ),
  ,
  drop = FALSE
]
estimated_offsets <- offset_histogram[
  !duplicated(offset_histogram$readlen),
  ,
  drop = FALSE
]
estimated_offsets$p_offset <- -estimated_offsets$distance
calibration_output <- Sys.getenv(
  "HEK293_OFFSET_CALIBRATION_OUTPUT",
  "ribostan_fixed_input_offset_calibration.tsv"
)
write.table(
  offset_histogram,
  calibration_output,
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)
if (identical(
  toupper(Sys.getenv("HEK293_OFFSET_CALIBRATION_ONLY", "FALSE")),
  "TRUE"
)) {
  print(estimated_offsets[c("readlen", "p_offset", "n")], row.names = FALSE)
  quit(save = "no", status = 0L)
}

add_occurrence <- function(records) {
  records <- records[order(
    records$read_name,
    records$transcript_id,
    records$position
  ), , drop = FALSE]
  key <- paste(records$read_name, records$transcript_id, sep = "\r")
  records$occurrence <- ave(
    seq_len(nrow(records)),
    key,
    FUN = seq_along
  )
  records
}

rt_records <- add_occurrence(tx_psite_records)
rs_records <- add_occurrence(rs_records)
record_delta <- merge(
  rt_records,
  rs_records,
  by = c("read_name", "transcript_id", "occurrence"),
  all = TRUE,
  suffixes = c("_ribotaper", "_ribostan")
)
record_delta$status <- ifelse(
  is.na(record_delta$position_ribotaper),
  "ribostan_only",
  ifelse(
    is.na(record_delta$position_ribostan),
    "ribotaper_only",
    ifelse(
      record_delta$position_ribotaper == record_delta$position_ribostan,
      "exact",
      "position_difference"
    )
  )
)
record_delta$position_difference <- record_delta$position_ribostan -
  record_delta$position_ribotaper
strand_by_tx <- vapply(
  exons_by_tx,
  function(exons) as.character(strand(exons)[1L]),
  character(1L)
)
transcript_strand <- strand_by_tx[match(
  record_delta$transcript_id,
  names(exons_by_tx)
)]
record_delta$transcript_strand <- transcript_strand
write.table(
  record_delta,
  record_output,
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

rs_positions_by_tx <- split(rs_records$position, rs_records$transcript_id)
build_from_positions <- function(orf_id, positions_index) {
  meta <- rt[match(orf_id, rt$orf_id_tr), , drop = FALSE]
  orf_start <- as.integer(meta$tx_start)
  orf_end <- as.integer(meta$tx_end)
  positions <- positions_index[[meta$transcript_id]]
  positions <- positions[positions >= orf_start & positions <= orf_end]
  as.integer(tabulate(
    positions - orf_start + 1L,
    nbins = orf_end - orf_start + 1L
  ))
}
rs_vectors <- setNames(
  lapply(orf_ids, build_from_positions, positions_index = rs_positions_by_tx),
  orf_ids
)

evaluate_vector <- function(vector) {
  stats <- run_reference(vector)
  n_psites <- sum(vector)
  frame <- if (sum(vector > 0L) <= 1L) {
    NA_real_
  } else {
    sum(vector[seq.int(1L, length(vector), by = 3L)]) / n_psites
  }
  call <- isTRUE(!is.na(stats["p.value"]) &
    stats["p.value"] < 0.05 &
    frame > 0.5 &
    n_psites >= 6L)
  c(
    n_psites = n_psites,
    p = unname(stats["p.value"]),
    frame = frame,
    call = call
  )
}
rs_stats <- t(vapply(rs_vectors, evaluate_vector, numeric(4L)))

orf_delta <- data.frame(
  orf_id_tr = orf_ids,
  vector_identical = vapply(
    orf_ids,
    function(orf_id) identical(vectors[[orf_id]], rs_vectors[[orf_id]]),
    logical(1L)
  ),
  n_ribotaper = result$n_psites,
  n_ribostan = rs_stats[, "n_psites"],
  p_ribotaper = result$p_ribotaper,
  p_ribostan = rs_stats[, "p"],
  frame_ribotaper = result$pct_inframe_ribotaper,
  frame_ribostan = rs_stats[, "frame"],
  call_ribotaper = result$call_ribotaper,
  call_ribostan = as.logical(rs_stats[, "call"]),
  stringsAsFactors = FALSE
)
orf_delta$outcome <- ifelse(
  orf_delta$call_ribotaper & !orf_delta$call_ribostan,
  "lost",
  ifelse(
    !orf_delta$call_ribotaper & orf_delta$call_ribostan,
    "gained",
    "unchanged"
  )
)
write.table(
  orf_delta,
  orf_output,
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

cat(
  if (nzchar(offset_path)) {
    "\nSTAGE 3: estimated-offset delta\n"
  } else {
    "\nSTAGE 2: P-site construction only\n"
  }
)
cat(
  "Offsets:",
  paste(
    sprintf(
      "%d:%d",
      sort(unique(selection_width)),
      vapply(
        sort(unique(selection_width)),
        function(read_length) unique(offset_by_record[selection_width == read_length]),
        integer(1L)
      )
    ),
    collapse = ", "
  ),
  "\n"
)
cat("Fixed-input metagene peak estimates:\n")
print(estimated_offsets[c("readlen", "p_offset", "n")], row.names = FALSE)
cat("Per-read/transcript P-site records:\n")
print(table(record_delta$status))
cat("\nRecords by transcript strand:\n")
print(table(record_delta$status, record_delta$transcript_strand))
shared <- !is.na(record_delta$position_ribotaper) &
  !is.na(record_delta$position_ribostan)
cat(sprintf(
  "Exact positions among shared records: %d/%d\n",
  sum(record_delta$status[shared] == "exact"),
  sum(shared)
))
cat("\nPer-ORF vectors:\n")
print(table(orf_delta$vector_identical))
cat("\nCall outcomes:\n")
print(table(orf_delta$outcome))
message("Wrote ", record_output, " and ", orf_output)
