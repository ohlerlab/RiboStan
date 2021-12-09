test_that("estimating per-codon occupanies works", 
{
    #first, let's load up some test data
  anno_file <- here::here('test.gc32.gtf')
  if(!file.exists(anno_file)){
    library(AnnotationHub)
    ah <- AnnotationHub()
    gencode32 <- ah[['AH75191']]
    seqlevels(gencode32)<-'chr22'
    suppressWarnings({rtracklayer::export(gencode32, anno_file, format='GTF')})
  }
  fafile <- here::here('chr22.fa')
  library(BSgenome.Hsapiens.UCSC.hg38)
  if(!file.exists(fafile)){
    seq <- Biostrings::DNAStringSet(BSgenome.Hsapiens.UCSC.hg38[['chr22']])
   names(seq)<-'chr22'
   Biostrings::writeXStringSet(
    seq, fafile)
  }
  # file.copy(fafile, ".")
  # system("gunzip -f chr22.fa.gz")
  # fafile <- "chr22.fa"
  chr22_anno <- load_annotation(anno_file, fafile)
  data(rpfs)
  data(offsets_df)
  psites <- get_psite_gr(rpfs, offsets_df, chr22_anno)
  rust_codon_occ_df <- get_codon_occs(psites, offsets_df, chr22_anno,
    n_genes = 500, method = "RUST"
  )
  codon_occ_dfcols <- c("codon", "position", "estimate", "upper", "lower", "p.value")
  expect_equal(colnames(rust_codon_occ_df), codon_occ_dfcols)
  expect_true(all(is.na(rust_codon_occ_df$upper)))
  expect_true(all(is.na(rust_codon_occ_df$lower)))
  expect_true(all(!is.na(rust_codon_occ_df$estimate)))
  codons <- Biostrings::GENETIC_CODE
  nonstops <- names(codons)[codons != "*"]
  expect_true(all(nonstops %in% rust_codon_occ_df$codon))
  expect_true(all(c("a_codon", "p_codon") %in% rust_codon_occ_df$position))
  glm_rust_codon_occ_df <- get_codon_occs(psites, offsets_df, chr22_anno,
    n_genes = 500, method = "RUST_glm"
  )
  tr_elong <- get_orf_elong(chr22_anno, glm_rust_codon_occ_df)
  expect_true(all(names(chr22_anno$cdsgrl) %in% tr_elong$orf_id))
  expect_equal(colnames(tr_elong), c("orf_id", "mean_occ"))
}
)
