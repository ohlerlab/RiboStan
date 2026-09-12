# Shared HEK293 benchmark annotation resolution (RiboTaper-style GTF filtering).
#
# RiboTaper's create_annotations_files.bash defaults for GENCODE v19 are
# use_ccdsid=true and use_appris=true on the full gencode.v19.annotation.gtf.
# Set HEK293_GTF_NAME to bypass auto-build and use a pre-built GTF path.

#' Parse one GTF attribute value from a line.
#' @keywords internal
.gtf_attr <- function(line, key) {
  m <- regmatches(
    line,
    regexec(paste0(key, ' "([^"]+)"'), line, perl = TRUE)
  )[[1L]]
  if (length(m) == 0L) {
    return(NA_character_)
  }
  m[[2L]]
}

#' Filter a GENCODE GTF to RiboTaper's CCDS/APPRIS transcript selection.
#'
#' Replicates the transcript-bucketing in RiboTaper
#' \code{create_annotations_files.bash} (ohlerlab/RiboTaper): CCDS-tagged CDS
#' rows define CCDS genes; within those genes APPRIS picks one isoform set
#' per gene (principal first, else any APPRIS, else all isoforms); non-CCDS
#' genes keep every transcript.
#'
#' @param source_gtf Path to the full source GTF (typically
#'   \code{gencode.v19.annotation.gtf}).
#' @param dest_gtf Output GTF path.
#' @param use_ccds If \code{TRUE}, CCDS genes are those with a \code{ccdsid}
#'   tag on a CDS row (RiboTaper \code{use_ccdsid=true}). If \code{FALSE}, any
#'   CDS-bearing gene counts as CCDS.
#' @param use_appris If \code{TRUE}, apply RiboTaper's one-isoform-per-CCDS-gene
#'   APPRIS rule. If \code{FALSE}, keep all transcripts in CCDS genes.
#' @param summary_path Optional path for a tab-separated summary of counts.
#' @return A list with kept transcript/gene counts and output paths.
#' @keywords internal
filter_gtf_ribotaper_style <- function(source_gtf,
                                       dest_gtf,
                                       use_ccds = TRUE,
                                       use_appris = TRUE,
                                       summary_path = NULL) {
  stopifnot(file.exists(source_gtf))
  lines <- readLines(source_gtf, warn = FALSE)
  header <- lines[grepl("^#", lines)]
  body <- lines[!grepl("^#", lines) & nzchar(trimws(lines))]
  if (length(body) == 0L) {
    stop("No data rows in GTF: ", source_gtf)
  }

  types <- vapply(strsplit(body, "\t", fixed = TRUE), `[`, character(1), 3)
  gene_ids <- vapply(body, .gtf_attr, character(1), key = "gene_id")
  tx_ids <- vapply(body, .gtf_attr, character(1), key = "transcript_id")

  cds_idx <- which(types == "CDS")
  if (length(cds_idx) == 0L) {
    stop("GTF has no CDS rows: ", source_gtf)
  }
  if (use_ccds) {
    genes_ccds <- unique(gene_ids[cds_idx[grepl("ccdsid", body[cds_idx], ignore.case = TRUE)]])
  } else {
    genes_ccds <- unique(gene_ids[cds_idx])
  }
  genes_ccds <- genes_ccds[!is.na(genes_ccds)]

  exon_idx <- which(types == "exon" & !is.na(tx_ids))
  tx_gene <- gene_ids[exon_idx]
  names(tx_gene) <- tx_ids[exon_idx]
  tx_gene <- tx_gene[!duplicated(names(tx_gene))]

  tx_lines <- split(body[exon_idx], tx_ids[exon_idx])
  tx_has_prin <- vapply(tx_lines, function(x) any(grepl("appris_prin", x)), logical(1))
  tx_has_appris <- vapply(tx_lines, function(x) any(grepl("appris", x, ignore.case = TRUE)), logical(1))

  all_genes <- unique(gene_ids[!is.na(gene_ids)])
  non_ccds_genes <- setdiff(all_genes, genes_ccds)

  kept_tx <- names(tx_gene)[tx_gene %in% non_ccds_genes]
  for (gid in genes_ccds) {
    txs <- names(tx_gene)[tx_gene == gid]
    if (length(txs) == 0L) {
      next
    }
    if (!use_appris) {
      kept_tx <- c(kept_tx, txs)
      next
    }
    prin <- txs[tx_has_prin[txs]]
    if (length(prin) > 0L) {
      kept_tx <- c(kept_tx, prin)
    } else {
      appr <- txs[tx_has_appris[txs]]
      if (length(appr) > 0L) {
        kept_tx <- c(kept_tx, appr)
      } else {
        kept_tx <- c(kept_tx, txs)
      }
    }
  }
  kept_tx <- unique(kept_tx)
  kept_genes <- unique(c(non_ccds_genes, genes_ccds))

  keep <- !is.na(gene_ids) & gene_ids %in% kept_genes &
    (is.na(tx_ids) | tx_ids %in% kept_tx)
  out_lines <- c(header, body[keep])
  writeLines(out_lines, dest_gtf)

  tx_whitelist <- file.path(dirname(dest_gtf), "ccds_appris_transcript_whitelist.txt")
  gene_whitelist <- file.path(dirname(dest_gtf), "ccds_appris_gene_whitelist.txt")
  writeLines(sort(kept_tx), tx_whitelist)
  writeLines(sort(kept_genes), gene_whitelist)

  n_tx_in <- length(unique(tx_ids[!is.na(tx_ids)]))
  summary <- data.frame(
    source_gtf = basename(source_gtf),
    dest_gtf = basename(dest_gtf),
    use_ccds = use_ccds,
    use_appris = use_appris,
    n_genes_ccds = length(genes_ccds),
    n_genes_non_ccds = length(non_ccds_genes),
    n_transcripts_in = n_tx_in,
    n_transcripts_kept = length(kept_tx),
    n_genes_kept = length(kept_genes),
    n_gtf_rows_out = length(body[keep]),
    stringsAsFactors = FALSE
  )
  if (!is.null(summary_path)) {
    write.table(summary, summary_path, sep = "\t", row.names = FALSE, quote = FALSE)
  }

  message(sprintf(
    "RiboTaper-style GTF filter: %s -> %s (%d/%d transcripts, use_ccds=%s use_appris=%s)",
    basename(source_gtf), basename(dest_gtf),
    length(kept_tx), n_tx_in, use_ccds, use_appris
  ))

  list(
    path = dest_gtf,
    source = source_gtf,
    use_ccds = use_ccds,
    use_appris = use_appris,
    summary = summary,
    transcript_whitelist = tx_whitelist,
    gene_whitelist = gene_whitelist
  )
}

#' Resolve the annotation GTF for HEK293 benchmark scripts.
#'
#' If \code{HEK293_GTF_NAME} is set, that file under \code{workdir} is used
#' as-is. Otherwise build (or reuse) a filtered GTF from \code{HEK293_GTF_SOURCE}
#' using \code{HEK293_USE_CCDS} / \code{HEK293_USE_APPRIS} (RiboTaper defaults:
#' both \code{TRUE} on GENCODE v19).
#'
#' @param workdir Benchmark working directory containing GTF/FASTA files.
#' @param rebuild If \code{TRUE}, regenerate the filtered GTF even when present.
#' @return A list with \code{path}, \code{name}, \code{source}, filter flags,
#'   and \code{built} (\code{TRUE} if filter ran this call).
resolve_hek293_gtf <- function(workdir, rebuild = FALSE) {
  explicit <- Sys.getenv("HEK293_GTF_NAME", "")
  if (nzchar(explicit)) {
    path <- file.path(workdir, explicit)
    return(list(
      path = path,
      name = explicit,
      source = NA_character_,
      use_ccds = NA,
      use_appris = NA,
      built = FALSE,
      explicit = TRUE
    ))
  }

  source_name <- Sys.getenv("HEK293_GTF_SOURCE", "gencode.v19.annotation.gtf")
  use_ccds <- as.logical(Sys.getenv("HEK293_USE_CCDS", "TRUE"))
  use_appris <- as.logical(Sys.getenv("HEK293_USE_APPRIS", "TRUE"))
  stopifnot(!is.na(use_ccds), !is.na(use_appris))

  dest_name <- if (use_ccds && use_appris) {
    "gencode.v19.ccds_appris.gtf"
  } else if (use_ccds) {
    "gencode.v19.ccds.gtf"
  } else if (use_appris) {
    "gencode.v19.appris.gtf"
  } else {
    basename(source_name)
  }

  source_path <- file.path(workdir, source_name)
  summary_path <- file.path(workdir, "hek293_gtf_filter_summary.txt")

  if (!use_ccds && !use_appris) {
    if (!file.exists(source_path)) {
      stop("Missing source GTF: ", source_path)
    }
    return(list(
      path = source_path,
      name = source_name,
      source = source_name,
      use_ccds = use_ccds,
      use_appris = use_appris,
      built = FALSE,
      explicit = FALSE,
      summary_path = summary_path
    ))
  }

  dest_path <- file.path(workdir, dest_name)
  if (rebuild || !file.exists(dest_path)) {
    if (!file.exists(source_path)) {
      stop("Missing source GTF for RiboTaper-style filter: ", source_path)
    }
    filter_gtf_ribotaper_style(
      source_path, dest_path,
      use_ccds = use_ccds,
      use_appris = use_appris,
      summary_path = summary_path
    )
    built <- TRUE
  } else {
    built <- FALSE
  }

  list(
    path = dest_path,
    name = dest_name,
    source = source_name,
    use_ccds = use_ccds,
    use_appris = use_appris,
    built = built,
    explicit = FALSE,
    summary_path = summary_path
  )
}
