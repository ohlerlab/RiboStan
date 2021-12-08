test_that("test optimisation of quantification works", 
{
  data(chr22_anno)
  data(rpfs)
  data(offsets_df)
  psites <- get_psite_gr(rpfs, offsets_df, chr22_anno)
  # now, use Stan to estimate normalized p-site densities for our data
  olduORFnum <- chr22_anno$uORF%>%sum
  filt_chr22_anno <- periodicity_filter_uORFs(psites, chr22_anno,
    remove=TRUE, n_cores=1)
  expect_true(sum(filt_chr22_anno$uORF) < sum(chr22_anno$uORF))
}
)
