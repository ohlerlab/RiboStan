# Unit tests for RiboTaper-style GTF filtering (HEK293 benchmark helper).

test_that("filter_gtf_ribotaper_style applies CCDS+APPRIS transcript selection per gene", {
  bench_dir <- normalizePath(
    file.path(testthat::test_path(), "..", "..", "scripts", "hek293_ribotaper_benchmark"),
    mustWork = TRUE
  )
  source(file.path(bench_dir, "hek293_gtf_config.R"), local = TRUE)

  tmp <- tempfile("gtf_filter_")
  dir.create(tmp)
  src <- file.path(tmp, "mini.gtf")
  out <- file.path(tmp, "mini.filtered.gtf")

  # Gene A (CCDS): principal + alternate isoforms -> keep principal only.
  # Gene B (CCDS): alternate APPRIS only -> keep alternate.
  # Gene C (CCDS): no APPRIS -> keep both isoforms.
  # Gene D (lncRNA, no CCDS): keep its transcript.
  lines <- c(
    paste("chr1", "Ens", "gene", "1", "1000", ".", "+", ".",
          'gene_id "G_A"; gene_type "protein_coding";', sep = "\t"),
    paste("chr1", "Ens", "transcript", "1", "1000", ".", "+", ".",
          'gene_id "G_A"; transcript_id "TX_A1"; tag "appris_principal";', sep = "\t"),
    paste("chr1", "Ens", "exon", "1", "1000", ".", "+", ".",
          'gene_id "G_A"; transcript_id "TX_A1"; tag "appris_principal";', sep = "\t"),
    paste("chr1", "Ens", "CDS", "10", "999", ".", "+", "0",
          'gene_id "G_A"; transcript_id "TX_A1"; tag "CCDS"; ccdsid "CCDS1";', sep = "\t"),
    paste("chr1", "Ens", "transcript", "1", "1000", ".", "+", ".",
          'gene_id "G_A"; transcript_id "TX_A2"; tag "appris_alternative_1";', sep = "\t"),
    paste("chr1", "Ens", "exon", "1", "1000", ".", "+", ".",
          'gene_id "G_A"; transcript_id "TX_A2"; tag "appris_alternative_1";', sep = "\t"),
    paste("chr1", "Ens", "CDS", "10", "999", ".", "+", "0",
          'gene_id "G_A"; transcript_id "TX_A2"; tag "CCDS"; ccdsid "CCDS1";', sep = "\t"),

    paste("chr1", "Ens", "gene", "1", "2000", ".", "+", ".",
          'gene_id "G_B"; gene_type "protein_coding";', sep = "\t"),
    paste("chr1", "Ens", "exon", "2000", "3000", ".", "+", ".",
          'gene_id "G_B"; transcript_id "TX_B1"; tag "appris_alternative_1";', sep = "\t"),
    paste("chr1", "Ens", "CDS", "2010", "2999", ".", "+", "0",
          'gene_id "G_B"; transcript_id "TX_B1"; ccdsid "CCDS2";', sep = "\t"),

    paste("chr1", "Ens", "gene", "1", "4000", ".", "+", ".",
          'gene_id "G_C"; gene_type "protein_coding";', sep = "\t"),
    paste("chr1", "Ens", "exon", "4000", "5000", ".", "+", ".",
          'gene_id "G_C"; transcript_id "TX_C1";', sep = "\t"),
    paste("chr1", "Ens", "exon", "4000", "5000", ".", "+", ".",
          'gene_id "G_C"; transcript_id "TX_C2";', sep = "\t"),
    paste("chr1", "Ens", "CDS", "4010", "4999", ".", "+", "0",
          'gene_id "G_C"; transcript_id "TX_C1"; ccdsid "CCDS3";', sep = "\t"),

    paste("chr1", "Ens", "gene", "6000", "7000", ".", "+", ".",
          'gene_id "G_D"; gene_type "lncRNA";', sep = "\t"),
    paste("chr1", "Ens", "exon", "6000", "7000", ".", "+", ".",
          'gene_id "G_D"; transcript_id "TX_D1";', sep = "\t")
  )
  writeLines(lines, src)

  res <- filter_gtf_ribotaper_style(src, out, use_ccds = TRUE, use_appris = TRUE)
  kept <- readLines(out, warn = FALSE)
  kept_tx <- unique(vapply(kept, .gtf_attr, character(1), key = "transcript_id"))
  kept_tx <- kept_tx[!is.na(kept_tx)]

  expect_setequal(kept_tx, c("TX_A1", "TX_B1", "TX_C1", "TX_C2", "TX_D1"))
  expect_equal(res$summary$n_transcripts_kept, 5L)
})
