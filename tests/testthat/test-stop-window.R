test_that("load_annotation retains the first uORF stop-codon nucleotide", {
  sequence <- rep("C", 100L)
  sequence[10:18] <- strsplit("ATGAAATAA", "", fixed = TRUE)[[1L]]
  sequence[50:58] <- strsplit("ATGAAATAA", "", fixed = TRUE)[[1L]]

  fasta <- tempfile(fileext = ".fa")
  dna <- Biostrings::DNAStringSet(paste(sequence, collapse = ""))
  names(dna) <- "chr1"
  Biostrings::writeXStringSet(dna, fasta)
  gtf <- tempfile(fileext = ".gtf")
  attributes <- 'gene_id "gene1"; transcript_id "tx1"; gene_name "GENE1";'
  writeLines(c(
    paste("chr1", "test", "gene", 1L, 100L, ".", "+", ".", attributes,
          sep = "\t"),
    paste("chr1", "test", "transcript", 1L, 100L, ".", "+", ".",
          attributes, sep = "\t"),
    paste("chr1", "test", "exon", 1L, 100L, ".", "+", ".", attributes,
          sep = "\t"),
    paste("chr1", "test", "CDS", 50L, 55L, ".", "+", "0", attributes,
          sep = "\t"),
    paste("chr1", "test", "start_codon", 50L, 52L, ".", "+", "0",
          attributes, sep = "\t"),
    paste("chr1", "test", "stop_codon", 56L, 58L, ".", "+", "0",
          attributes, sep = "\t")
  ), gtf)

  anno <- suppressWarnings(suppressMessages(load_annotation(gtf, fasta)))
  uorf_ids <- names(anno$uORF)[anno$uORF]

  expect_length(uorf_ids, 1L)
  expect_equal(as.integer(width(anno$trspacecds[uorf_ids])), 7L)
  expect_equal(as.integer(start(anno$trspacecds[uorf_ids])), 10L)
  expect_equal(as.integer(end(anno$trspacecds[uorf_ids])), 16L)
})

test_that("invalid main CDS does not suppress valid uORF discovery", {
  sequence <- rep("C", 300L)
  sequence[50:58] <- strsplit("ATGAAATAA", "", fixed = TRUE)[[1L]]
  sequence[150:158] <- strsplit("ATGAAATAA", "", fixed = TRUE)[[1L]]
  sequence[210:218] <- strsplit("ATGAAATAA", "", fixed = TRUE)[[1L]]
  sequence[250:258] <- strsplit("ATGAAAAAA", "", fixed = TRUE)[[1L]]

  fasta <- tempfile(fileext = ".fa")
  dna <- Biostrings::DNAStringSet(paste(sequence, collapse = ""))
  names(dna) <- "chr1"
  Biostrings::writeXStringSet(dna, fasta)
  gtf <- tempfile(fileext = ".gtf")
  transcript_lines <- function(transcript_id, gene_id, offset) {
    attributes <- sprintf(
      'gene_id "%s"; transcript_id "%s"; gene_name "%s";',
      gene_id, transcript_id, gene_id
    )
    paste(
      "chr1", "test",
      c("gene", "transcript", "exon", "CDS", "start_codon", "stop_codon"),
      offset + c(1L, 1L, 1L, 50L, 50L, 56L),
      offset + c(100L, 100L, 100L, 55L, 52L, 58L),
      ".", "+", c(".", ".", ".", "0", "0", "0"), attributes,
      sep = "\t"
    )
  }
  writeLines(c(
    transcript_lines("tx_valid1", "gene1", 0L),
    transcript_lines("tx_valid2", "gene2", 100L),
    transcript_lines("tx_invalid", "gene3", 200L)
  ), gtf)

  anno <- suppressWarnings(suppressMessages(load_annotation(gtf, fasta)))
  uorf_ids <- names(anno$uORF)[anno$uORF]
  invalid_uorf_ids <- uorf_ids[startsWith(uorf_ids, "tx_invalid_")]

  expect_length(invalid_uorf_ids, 1L)
  expect_false("tx_invalid" %in% names(anno$trspacecds))
  expect_equal(
    as.integer(start(anno$trspacecds[invalid_uorf_ids])),
    10L
  )
})
