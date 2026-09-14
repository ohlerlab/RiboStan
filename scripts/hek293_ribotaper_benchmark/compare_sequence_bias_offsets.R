#!/usr/bin/env Rscript
# Test whether footprint-end sequence predicts residual P-site phase errors.
args <- commandArgs(trailingOnly = TRUE)
mode <- if (length(args)) args[[1L]] else "synthetic"
workdir <- if (length(args) >= 2L) args[[2L]] else "."
stopifnot(mode %in% c("synthetic", "real", "both"))
set.seed(104L)

context_radius <- as.integer(Sys.getenv("HEK293_SEQUENCE_CONTEXT_RADIUS", "3"))
n_folds <- as.integer(Sys.getenv("HEK293_SEQUENCE_FOLDS", "5"))
max_train <- as.integer(Sys.getenv("HEK293_SEQUENCE_MAX_TRAIN", "200000"))
max_calibration <- as.integer(Sys.getenv(
  "HEK293_SEQUENCE_MAX_CALIBRATION",
  "500000"
))

phase_score <- function(frame) {
  counts <- tabulate(as.integer(frame) + 1L, nbins = 3L)
  sum((counts - mean(counts))^2) / sum(counts)^2
}

phase0_fraction <- function(frame) mean(as.integer(frame) == 0L)

fit_phase_model <- function(data, feature_cols) {
  formula <- stats::reformulate(feature_cols, response = "target")
  lapply(0:2, function(value) {
    model_data <- data
    model_data$target <- as.integer(model_data$correction) - 1L == value
    if (length(unique(model_data$target)) == 1L) {
      return(list(constant = mean(model_data$target)))
    }
    stats::glm(formula, data = model_data, family = stats::binomial())
  })
}

predict_phase_model <- function(model, data) {
  probabilities <- vapply(model, function(component) {
    if (is.list(component) && !is.null(component$constant)) {
      return(rep(component$constant, nrow(data)))
    }
    stats::predict(component, newdata = data, type = "response")
  }, numeric(nrow(data)))
  probabilities[!is.finite(probabilities)] <- 0
  prediction <- max.col(probabilities, ties.method = "first") - 1L
  prediction[rowSums(probabilities) == 0] <- 0L
  prediction
}

model_metrics <- function(data, prediction, label) {
  corrected_frame <- (data$observed_frame + prediction) %% 3L
  data.frame(
    model = label,
    accuracy = mean(prediction == as.integer(as.character(data$correction))),
    phase0_fraction = phase0_fraction(corrected_frame),
    phase_score = phase_score(corrected_frame),
    correction_rate = mean(prediction != 0L),
    stringsAsFactors = FALSE
  )
}

random_sequence <- function(n) {
  paste0(sample(c("A", "C", "G", "T"), n, replace = TRUE), collapse = "")
}

sequence_bases <- function(sequence, positions) {
  valid <- positions >= 1L & positions <= nchar(sequence)
  bases <- rep("N", length(positions))
  bases[valid] <- substring(sequence, positions[valid], positions[valid])
  bases
}

offset_label <- function(prefix, offset) {
  paste0(prefix, "_", ifelse(offset < 0L, "m", "p"), abs(offset))
}

add_synthetic_context <- function(data, sequences) {
  offsets <- -context_radius:context_radius
  for (offset in offsets) {
    data[[offset_label("five", offset)]] <- vapply(
      seq_len(nrow(data)),
      function(index) sequence_bases(
        sequences[[data$gene[[index]]]],
        data$read_start[[index]] + offset
      ),
      character(1L)
    )
    data[[offset_label("three", offset)]] <- vapply(
      seq_len(nrow(data)),
      function(index) sequence_bases(
        sequences[[data$gene[[index]]]],
        data$read_end[[index]] + offset
      ),
      character(1L)
    )
  }
  context_cols <- grep("^(five|three)_", names(data), value = TRUE)
  data[context_cols] <- lapply(
    data[context_cols],
    factor,
    levels = c("A", "C", "G", "T", "N")
  )
  data$readlen <- factor(data$readlen, levels = 26:30)
  data
}

simulate_cleavage_data <- function(n_genes = 240L, reads_per_gene = 120L,
                                   sequence_bias = TRUE) {
  genes <- sprintf("gene_%03d", seq_len(n_genes))
  sequences <- stats::setNames(
    replicate(n_genes, random_sequence(1200L)),
    genes
  )
  rows <- lapply(genes, function(gene) {
    true_psite <- sample(seq.int(100L, 1100L, by = 3L), reads_per_gene,
                         replace = TRUE)
    readlen <- sample(26:30, reads_per_gene, replace = TRUE)
    nominal_start <- true_psite - 12L
    motif <- sequence_bases(sequences[[gene]], nominal_start)
    displacement <- if (sequence_bias) {
      ifelse(motif == "G", 1L, ifelse(motif == "A", -1L, 0L))
    } else {
      rep(0L, reads_per_gene)
    }
    read_start <- nominal_start + displacement
    data.frame(
      gene = gene,
      readlen = readlen,
      read_start = read_start,
      read_end = read_start + readlen - 1L,
      observed_frame = displacement %% 3L,
      correction = factor((-displacement) %% 3L, levels = 0:2)
    )
  })
  add_synthetic_context(do.call(rbind, rows), sequences)
}

run_synthetic <- function(output_dir) {
  biased <- simulate_cleavage_data(sequence_bias = TRUE)
  unbiased <- simulate_cleavage_data(sequence_bias = FALSE)
  genes <- unique(biased$gene)
  test_genes <- sample(genes, ceiling(length(genes) * 0.2))
  train <- biased[!biased$gene %in% test_genes, ]
  test <- biased[biased$gene %in% test_genes, ]
  context_cols <- grep("^(five|three)_", names(train), value = TRUE)
  sequence_features <- c("readlen", context_cols)

  length_model <- fit_phase_model(train, "readlen")
  sequence_model <- fit_phase_model(train, sequence_features)
  permuted <- train
  permuted$correction <- sample(permuted$correction)
  permutation_model <- fit_phase_model(permuted, sequence_features)

  predictions <- list(
    baseline = rep(0L, nrow(test)),
    read_length = predict_phase_model(length_model, test),
    sequence = predict_phase_model(sequence_model, test),
    permuted_labels = predict_phase_model(permutation_model, test),
    oracle = as.integer(as.character(test$correction))
  )
  metrics <- do.call(rbind, Map(
    function(prediction, label) model_metrics(test, prediction, label),
    predictions,
    names(predictions)
  ))
  metrics$dataset <- "known_sequence_bias"

  unbiased_train <- unbiased[!unbiased$gene %in% test_genes, ]
  unbiased_test <- unbiased[unbiased$gene %in% test_genes, ]
  unbiased_model <- fit_phase_model(unbiased_train, sequence_features)
  no_bias_metrics <- model_metrics(
    unbiased_test,
    predict_phase_model(unbiased_model, unbiased_test),
    "sequence"
  )
  no_bias_metrics$dataset <- "no_sequence_bias"
  metrics <- rbind(metrics, no_bias_metrics)

  null_prediction <- predict_phase_model(sequence_model, test)
  baseline_null_fpr <- replicate(500L, {
    frame <- sample(0:2, 60L, replace = TRUE)
    stats::chisq.test(tabulate(frame + 1L, nbins = 3L))$p.value < 0.05
  })
  corrected_null_fpr <- replicate(500L, {
    index <- sample(seq_along(null_prediction), 60L)
    frame <- sample(0:2, 60L, replace = TRUE)
    corrected <- (frame + null_prediction[index]) %% 3L
    stats::chisq.test(tabulate(corrected + 1L, nbins = 3L))$p.value < 0.05
  })
  null_summary <- data.frame(
    model = c("baseline", "sequence"),
    metric = "uniform_phase_chisq_fpr",
    value = c(mean(baseline_null_fpr), mean(corrected_null_fpr))
  )

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  write.table(
    metrics,
    file.path(output_dir, "sequence_bias_synthetic_metrics.tsv"),
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )
  write.table(
    null_summary,
    file.path(output_dir, "sequence_bias_synthetic_null.tsv"),
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )
  print(metrics)
  print(null_summary)
  stopifnot(
    metrics$accuracy[metrics$dataset == "known_sequence_bias" &
      metrics$model == "sequence"] > 0.65,
    metrics$phase0_fraction[metrics$dataset == "known_sequence_bias" &
      metrics$model == "sequence"] >
      metrics$phase0_fraction[metrics$dataset == "known_sequence_bias" &
        metrics$model == "read_length"]
  )
}

extract_context <- function(data, transcript_sequences) {
  offsets <- -context_radius:context_radius
  context <- matrix(
    "N",
    nrow = nrow(data),
    ncol = length(offsets) * 2L
  )
  colnames(context) <- c(
    offset_label("five", offsets),
    offset_label("three", offsets)
  )
  by_transcript <- split(seq_len(nrow(data)), data$transcript_id)
  for (transcript_id in names(by_transcript)) {
    index <- by_transcript[[transcript_id]]
    sequence <- transcript_sequences[[transcript_id]]
    if (is.null(sequence)) {
      next
    }
    sequence <- as.character(sequence)
    for (offset_index in seq_along(offsets)) {
      offset <- offsets[[offset_index]]
      context[index, offset_index] <- sequence_bases(
        sequence,
        data$read_start[index] + offset
      )
      context[index, length(offsets) + offset_index] <- sequence_bases(
        sequence,
        data$read_end[index] + offset
      )
    }
  }
  context <- as.data.frame(context, stringsAsFactors = FALSE)
  context[] <- lapply(
    context,
    factor,
    levels = c("A", "C", "G", "T", "N")
  )
  data$readlen <- factor(data$readlen, levels = 26:30)
  cbind(data, context)
}

make_gene_folds <- function(genes) {
  unique_genes <- sample(unique(genes))
  stats::setNames((seq_along(unique_genes) - 1L) %% n_folds + 1L,
                  unique_genes)
}

sample_training_rows <- function(index) {
  if (length(index) > max_train) sample(index, max_train) else index
}

predict_crossfit <- function(models, data, folds) {
  prediction <- integer(nrow(data))
  for (fold in seq_along(models)) {
    index <- which(folds == fold)
    if (length(index)) {
      prediction[index] <- predict_phase_model(models[[fold]], data[index, ])
    }
  }
  prediction
}

call_table <- function(psites, anno, n_cores) {
  tested <- periodicity_filter_uORFs(
    psites,
    anno,
    remove = FALSE,
    min_psites = 6L,
    min_frame_frac = 0.5,
    n_cores = n_cores
  )
  ids <- names(anno$uORF)[anno$uORF]
  metadata <- as.data.frame(S4Vectors::mcols(tested$trspacecds[ids]))
  representative <- !is.na(metadata$locus_representative) &
    metadata$locus_representative
  data.frame(
    orf_id = ids[representative],
    orf_locus_id = metadata$orf_locus_id[representative],
    p.value = metadata$p.value[representative],
    q.value = metadata$q.value[representative],
    n_psites = metadata$n_psites[representative],
    pct_inframe = metadata$pct_inframe[representative],
    call = !is.na(metadata$q.value[representative]) &
      metadata$q.value[representative] < 0.05 &
      metadata$n_psites[representative] >= 6L &
      metadata$pct_inframe[representative] > 0.5,
    stringsAsFactors = FALSE
  )
}

run_real <- function(workdir, output_dir) {
  suppressPackageStartupMessages({
    library(Ribostan)
    library(GenomicRanges)
    library(GenomicFeatures)
  })
  script_dir <- dirname(sub(
    "--file=", "",
    commandArgs(trailingOnly = FALSE)[grep(
      "--file=", commandArgs(trailingOnly = FALSE)
    )]
  ))
  source(file.path(script_dir, "hek293_gtf_config.R"), local = TRUE)
  old_workdir <- setwd(workdir)
  on.exit(setwd(old_workdir), add = TRUE)

  gtf <- resolve_hek293_gtf(workdir)$path
  fasta <- file.path(workdir, "GRCh37.p13.genome.fa")
  projected_rds <- Sys.getenv(
    "HEK293_PROJECTED_RDS",
    "rpfs_genome_projected_rrnafilt_matched.rds"
  )
  offset_path <- Sys.getenv(
    "HEK293_SEQUENCE_OFFSETS",
    "ribostan_offsetsribotaper_matched.tsv"
  )
  n_cores <- as.integer(Sys.getenv("HEK293_N_CORES", "8"))

  anno_cds <- load_annotation(gtf, fasta, add_uorfs = FALSE)
  rpfs <- readRDS(projected_rds)
  offsets <- if (file.exists(offset_path)) {
    read.delim(offset_path, stringsAsFactors = FALSE)
  } else {
    calculated <- get_offsets(rpfs, anno_cds, method = "metagene_peak")
    write.table(
      calculated,
      offset_path,
      sep = "\t",
      row.names = FALSE,
      quote = FALSE
    )
    calculated
  }
  cds <- anno_cds$trspacecds
  hits <- findOverlaps(rpfs, cds, ignore.strand = TRUE)
  first_hit <- !duplicated(queryHits(hits))
  read_index <- queryHits(hits)[first_hit]
  cds_index <- subjectHits(hits)[first_hit]
  phase <- (start(rpfs)[read_index] - start(cds)[cds_index]) %% 3L
  lookup <- match(
    paste(rpfs$readlen[read_index], phase),
    paste(offsets$readlen, offsets$phase)
  )
  keep <- !is.na(lookup)
  read_index <- read_index[keep]
  cds_index <- cds_index[keep]
  phase <- phase[keep]
  p_offset <- offsets$p_offset[lookup[keep]]
  transcript_id <- as.character(seqnames(rpfs)[read_index])
  tx_gene <- anno_cds$trgiddf$gene_id
  names(tx_gene) <- anno_cds$trgiddf$orf_id
  gene <- unname(tx_gene[transcript_id])
  gene[is.na(gene)] <- transcript_id[is.na(gene)]
  calibration <- data.frame(
    read_index = read_index,
    transcript_id = transcript_id,
    gene = gene,
    readlen = as.integer(rpfs$readlen[read_index]),
    read_start = as.integer(start(rpfs)[read_index]),
    read_end = as.integer(start(rpfs)[read_index] +
      rpfs$readlen[read_index] - 1L),
    observed_frame = (phase + p_offset) %% 3L,
    correction = factor((-(phase + p_offset)) %% 3L, levels = 0:2)
  )
  if (nrow(calibration) > max_calibration) {
    calibration <- calibration[sample(nrow(calibration), max_calibration), ]
  }
  transcript_ids <- intersect(
    unique(calibration$transcript_id),
    names(anno_cds$exonsgrl)
  )
  transcript_sequences <- extractTranscriptSeqs(
    anno_cds$fafileob,
    anno_cds$exonsgrl[transcript_ids]
  )
  calibration <- extract_context(calibration, transcript_sequences)
  context_cols <- grep("^(five|three)_", names(calibration), value = TRUE)
  sequence_features <- c("readlen", context_cols)
  gene_folds <- make_gene_folds(calibration$gene)
  calibration$fold <- unname(gene_folds[calibration$gene])
  calibration$pred_length <- NA_integer_
  calibration$pred_sequence <- NA_integer_
  calibration$pred_permuted <- NA_integer_
  sequence_models <- vector("list", n_folds)
  length_models <- vector("list", n_folds)
  permutation_models <- vector("list", n_folds)

  for (fold in seq_len(n_folds)) {
    train_index <- which(calibration$fold != fold)
    train_index <- sample_training_rows(train_index)
    test_index <- which(calibration$fold == fold)
    train <- calibration[train_index, ]
    sequence_models[[fold]] <- fit_phase_model(train, sequence_features)
    length_models[[fold]] <- fit_phase_model(train, "readlen")
    permuted <- train
    permuted$correction <- sample(permuted$correction)
    permutation_models[[fold]] <- fit_phase_model(
      permuted, sequence_features
    )
    calibration$pred_length[test_index] <- predict_phase_model(
      length_models[[fold]], calibration[test_index, ]
    )
    calibration$pred_sequence[test_index] <- predict_phase_model(
      sequence_models[[fold]], calibration[test_index, ]
    )
    calibration$pred_permuted[test_index] <- predict_phase_model(
      permutation_models[[fold]], calibration[test_index, ]
    )
  }

  real_metrics <- rbind(
    model_metrics(calibration, rep(0L, nrow(calibration)), "baseline"),
    model_metrics(calibration, calibration$pred_length, "read_length"),
    model_metrics(calibration, calibration$pred_sequence, "sequence"),
    model_metrics(calibration, calibration$pred_permuted, "permuted_labels")
  )
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  write.table(
    real_metrics,
    file.path(output_dir, "sequence_bias_real_metrics.tsv"),
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )

  reference <- read.delim(
    Sys.getenv("HEK293_PARITY_UORFS", "ribotaper_hek293_uorfs.tsv"),
    stringsAsFactors = FALSE
  )
  imported_gtf <- rtracklayer::import(gtf, format = "gtf")
  target_genes <- unique(as.character(
    imported_gtf$gene_id[
      as.character(imported_gtf$transcript_id) %in%
        unique(reference$transcript_id)
    ]
  ))
  keep_gtf <- as.character(imported_gtf$gene_id) %in% target_genes
  target_gtf <- tempfile(fileext = ".gtf")
  rtracklayer::export(imported_gtf[keep_gtf], target_gtf, format = "gtf")
  anno_full <- load_annotation(target_gtf, fasta, add_uorfs = TRUE)
  baseline_psites <- get_psite_gr(
    rpfs, offsets, anno_full, apply_phase_shiftdf = FALSE
  )
  psite_transcript <- as.character(seqnames(baseline_psites))
  psite_gene <- unname(tx_gene[psite_transcript])
  psite_gene[is.na(psite_gene)] <- psite_transcript[is.na(psite_gene)]
  psite_data <- data.frame(
    transcript_id = psite_transcript,
    gene = psite_gene,
    readlen = as.integer(baseline_psites$readlen),
    read_start = as.integer(start(baseline_psites) -
      baseline_psites$p_offset),
    read_end = as.integer(start(baseline_psites) -
      baseline_psites$p_offset + baseline_psites$readlen - 1L)
  )
  psite_data <- extract_context(psite_data, transcript_sequences)
  psite_fold <- unname(gene_folds[psite_data$gene])
  corrections <- list(
    read_length = predict_crossfit(length_models, psite_data, psite_fold),
    sequence = predict_crossfit(sequence_models, psite_data, psite_fold),
    permuted_labels = predict_crossfit(
      permutation_models, psite_data, psite_fold
    )
  )
  corrected_psites <- lapply(corrections, function(correction) {
    result <- shift(baseline_psites, correction)
    result[!Ribostan:::is_out_of_bounds(result)]
  })

  baseline_calls <- call_table(baseline_psites, anno_full, n_cores)
  corrected_calls <- lapply(
    corrected_psites,
    call_table,
    anno = anno_full,
    n_cores = n_cores
  )
  call_delta <- merge(
    baseline_calls,
    corrected_calls$sequence,
    by = "orf_locus_id",
    suffixes = c("_baseline", "_sequence"),
    all = TRUE
  )
  call_delta$call_baseline[is.na(call_delta$call_baseline)] <- FALSE
  call_delta$call_sequence[is.na(call_delta$call_sequence)] <- FALSE
  write.table(
    call_delta,
    file.path(output_dir, "sequence_bias_real_call_delta.tsv"),
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )
  call_sets <- c(list(baseline = baseline_calls), corrected_calls)
  for (label in names(call_sets)) {
    write.table(
      call_sets[[label]],
      file.path(output_dir, paste0("sequence_bias_real_calls_", label, ".tsv")),
      sep = "\t",
      row.names = FALSE,
      quote = FALSE
    )
  }
  call_summary <- do.call(rbind, lapply(names(call_sets), function(label) {
    data.frame(
      model = label,
      calls = sum(call_sets[[label]]$call),
      stringsAsFactors = FALSE
    )
  }))
  start_call_path <- Sys.getenv(
    "HEK293_SEQUENCE_START_CALLS",
    "ribostan_start_codon_delta.tsv"
  )
  if (file.exists(start_call_path)) {
    start_calls <- read.delim(start_call_path, stringsAsFactors = FALSE)
    atg <- start_calls[
      start_calls$scenario == "ATG" & start_calls$start_codon == "ATG",
    ]
    published_key <- paste(
      reference$transcript_id, reference$tx_start, reference$tx_end
    )
    atg_key <- paste(atg$transcript_id, atg$tx_start, atg$tx_end)
    published_loci <- unique(atg$orf_locus_id[atg_key %in% published_key])
    call_summary$published <- vapply(
      call_sets,
      function(result) sum(result$call &
        result$orf_locus_id %in% published_loci),
      integer(1L)
    )
    call_summary$nonpublished <- call_summary$calls -
      call_summary$published
  }
  write.table(
    call_summary,
    file.path(output_dir, "sequence_bias_real_call_summary.tsv"),
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )
  print(real_metrics)
  print(call_summary)
  print(with(call_delta, table(call_baseline, call_sequence)))
}

output_dir <- Sys.getenv(
  "HEK293_SEQUENCE_OUTPUT_DIR",
  file.path(workdir, "sequence_bias_results")
)
if (mode %in% c("synthetic", "both")) {
  run_synthetic(output_dir)
}
if (mode %in% c("real", "both")) {
  run_real(workdir, output_dir)
}
