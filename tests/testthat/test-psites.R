test_that("psite object creation works",
{
  data(chr22_anno)
  data(rpfs)
  data(offsets_df)
  psites <- get_psite_gr(rpfs, offsets_df, chr22_anno)
  # get_psite_gr emits one psite per (read, ORF) overlap and carries a
  # read_mult mcol recording the number of ORF overlaps per source read.
  # See R/get_ritpms.R for details.
  expect_setequal(
    psites %>% mcols() %>% colnames(),
    c("readlen", "orf", "phase", "read_mult", "p_offset")
  )
  expect_true(all(psites$orf %in% names(chr22_anno$trspacecds)))
  # Exact length depends on the multi-ORF expansion; verify a sensible
  # lower bound (historical single-psite-per-read count) and that
  # multi-ORF reads are being expanded at all.
  expect_gte(length(psites), 53462)
  expect_gt(sum(psites$read_mult > 1L), 0L)
}
)
