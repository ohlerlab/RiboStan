test_that("estimating per-codon occupanies works", 
{
  #first, let's load up some test data
  anno_file <- 'test.gc32.gtf'
    library(AnnotationHub)
    ah <- AnnotationHub()
    gencode32 <- ah[['AH75191']]
    seqlevels(gencode32)<-'chr22'
    rtracklayer::export(gencode32, anno_file, format='GTF')
  fafile <- 'chr22.fa'
  library(BSgenome.Hsapiens.UCSC.hg38)
    seq <- Biostrings::DNAStringSet(BSgenome.Hsapiens.UCSC.hg38[['chr22']])
   names(seq) <- 'chr22'
   Biostrings::writeXStringSet(
    seq, fafile)
  chr22_anno <- load_annotation(anno_file, fafile)
  data(rpfs)
  data(offsets_df)
  covgrs <- list(sample1 = rpfs)
  metacodondf <- get_metacodon_profs(covgrs, chr22_anno)
  expect_equal(colnames(metacodondf), c(
    "sample", "readlen", "codon",
    "position", "ro_cl", "re_c",
    "count", "nreadlen"
  ))
  rm(metacodondf)
  data(metacodondf)
  expect_equal(metacodondf$codon %>% table() %>% length(), 61)
  expect_equal(metacodondf %>% nrow(), 23058)
  expect_equal(metacodondf$ro_cl %>% is.na() %>% sum(), 0)
  expect_equal(metacodondf$re_c %>% is.na() %>% sum(), 0)
  expect_equal(metacodondf$count %>% is.na() %>% sum(), 0)
  # note this doesn't work that well on a small subset
  kl_df <- get_kl_df(metacodondf, chr22_anno)
  expect_true(!all(!is.finite(kl_df$KL)))
  expect_equal(
    colnames(kl_df),
    c("sample", "nreadlen", "position", "KL")
  )
  expect_equal(kl_df$nreadlen %>% table() %>% as.numeric() %>% unique(), 54)
  kl_offsets <- select_offsets(kl_df)
  expect_equal(colnames(kl_offsets), c("sample", "nreadlen", "p_offset"))
  expect_equal(nrow(kl_offsets), 7)
  allcodondt <- export_codon_dts(metacodondf, kl_offsets)
  expect_equal(
    allcodondt %>% colnames(),
    c("sample", "codon", "a_p3_site", "a_site", "e_site", "p_site")
  )
  expect_equal(allcodondt %>% nrow(), 61)
}
)
