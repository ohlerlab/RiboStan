test_that("test optimisation of quantification works", {
  data(chr22_anno)
  data(rpfs)
  data(offsets_df)
  data(ms_df)
  
  psites <- get_psite_gr(rpfs, offsets_df, chr22_anno)
  #now, use Stan to estimate normalized p-site densities for our data
  ritpms <- get_ritpms(psites, chr22_anno)
  #and get these at the gene level (ignoring uORFs)
  gritpms = gene_level_expr(ritpms, chr22_anno)
  #compare these to mass spec data
  gritpms <- gritpms%>%mutate_at('gene_id',str_replace,'\\.\\d+$','')
  compdf <- ms_df%>%left_join(gritpms)%>%
    filter(is.finite(log2(expr)))%>%
    filter(ribo>0)
  expect_true(nrow(compdf)==81)
  expect_true(cor(compdf$ribo,log2(compdf$expr))>0.75)


  ftests = ftest_orfs(psites%>%head(10000), chr22_anno) 
  expect_equal(colnames(ftests),c("orf_id", "spec_coef", "p.value"))
  
})