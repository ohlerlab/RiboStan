test_that("psite object creation works", {
  data(chr22_anno)
  data(rpfs)
  data(offsets_df)
  psites <- get_psite_gr(rpfs, offsets_df, chr22_anno)
  expect_equal(
    psites %>% mcols() %>% colnames(),
    c("readlen", "orf", "phase", "p_offset")
  )
  expect_true(all(psites$orf %in% names(chr22_anno$trspacecds)))
  expect_equal(length(psites), 63206)
})
