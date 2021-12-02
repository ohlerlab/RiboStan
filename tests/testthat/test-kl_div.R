test_that("estimating per-codon occupanies works", {
  data(chr22_anno)
  data(rpfs)
  data(offsets_df)
  covgrs = list(sample1=rpfs)
  metacodondf <- get_metacodon_profs(covgrs, chr22_anno) 
  expect_equal(colnames(metacodondf),c("sample", "readlen", "codon",
   "position", "ro_cl", "re_c",
  "count", "nreadlen"))
  expect_equal(metacodondf$codon%>%table%>%length,61)
  expect_equal(metacodondf%>%nrow,23058)
  expect_equal(metacodondf$ro_cl%>%is.na%>%sum,0)
  expect_equal(metacodondf$re_c%>%is.na%>%sum,0)
  expect_equal(metacodondf$count%>%is.na%>%sum,0)
  #note this doesn't work that well on a small subset
  kl_df<-get_kl_df(metacodondf, chr22_anno)
  expect_true(!all(!is.finite(kl_df$KL)))
  expect_equal(colnames(kl_df),
    c("sample", "nreadlen", "position", "KL"))
  expect_equal(kl_df$nreadlen%>%table%>%as.numeric%>%unique,54)
  kl_offsets <- select_offsets(kl_df)
  expect_equal(colnames(kl_offsets),c("sample", "nreadlen", "p_offset"))
  expect_equal(nrow(kl_offsets),5)
  allcodondt <- export_codon_dts(metacodondf, kl_offsets)
  expect_equal(allcodondt%>%colnames,
    c("sample", "codon", "a_p3_site", "a_site", "e_site", "p_site"))
  expect_equal(allcodondt%>%nrow,61)
 
})

