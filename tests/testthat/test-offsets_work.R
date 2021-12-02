test_that("finding offsets works", {
  data(chr22_anno)
  data(rpfs)
  data(offsets_df)
  testoffsets_df <- get_offsets(rpfs, chr22_anno)
  expect_identical(offsets_df, testoffsets_df)
})
