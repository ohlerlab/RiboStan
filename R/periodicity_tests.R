#' @importFrom multitaper dpss spec.mtm dropFreqs
NULL

#' Test a numeric vector for periodicity using the multitaper package
#'
#' @keywords Ribostan
#' @author Dermot Harnett, \email{dermot.p.harnett@gmail.com}
#'
#' @param psit A list of SimpleRleLists - RPF coverage split by readlength
#' @param k - The number of slepians to use for the test
#' @param bw - The bandwidth to use for the test
#'
#' @details This function pads short vectors up to 50, and performs a multitaper test for 3bp
#' peridicity on it's input, returning a spectral coefficient and p-value
#' @return a numeric vector with the spectral coefficient at 0.333... and the pvalue for the test

ftestvect <- function(psit, k = 24, bw = 12) {
  psit <- as.numeric(psit)
  slepians_values <- dpss(n = length(psit) %>% ifelse(. < 25, 50, .), k = k, nw = bw)
  # vals<-take_Fvals_spect(x = psit,n_tapers = k,time_bw = bw,slepians_values = sl)

  if (length(psit) < 25) {
    remain <- 50 - length(psit)
    halfrmn <- as.integer(remain / 2)
    psit <- c(rep(0, halfrmn), psit, rep(0, remain %% 2 + halfrmn))
  }
  #
  if (length(psit) < 1024 / 2) {
    padding <- 1024
  }
  if (length(psit) >= 1024 / 2) {
    padding <- "default"
  }
  #
  resSpec1 <- spec.mtm(as.ts(psit),
    k = k, nw = bw, nFFT = padding,
    centreWithSlepians = TRUE, Ftest = TRUE,
    maxAdaptiveIterations = 100, returnZeroFreq = FALSE,
    plot = FALSE, dpssIN = slepians_values
  )
  psit
  #
  resSpec2 <- dropFreqs(resSpec1, 0.29, 0.39)
  #
  closestfreqind <- which(abs((resSpec1$freq - (1 / 3))) == min(abs((resSpec1$freq - (1 / 3)))))
  #
  freq_max_3nt <- resSpec1$freq[closestfreqind]
  Fmax_3nt <- resSpec1$mtm$Ftest[closestfreqind]
  spect_3nt <- resSpec1$spec[closestfreqind]
  return(c(Fmax_3nt, spect_3nt))

  pval <- pf(q = vals[1], df1 = 2, df2 = (2 * 24) - 2, lower.tail = FALSE)
  return(c(vals[2], pval))
}


#' Run Multitaper Tests on All ORFs
#'
#' @keywords Ribostan
#' @author Dermot Harnett, \email{dermot.p.harnett@gmail.com}
#'
#' @param psites A GRanges object containing psites
#' @param anno annotation object
#' @param n_cores number of cores to use
#'
#' @details This function applies a multitaper test to
#' @return a numeric vector with the spectral coefficient at 0.333...
#' and the pvalue for the test
#' @export
#' @examples
#' data(chr22_anno)
#' data(rpfs)
#' data(offsets_df)
#' psites <- get_psite_gr(rpfs, offsets_df, chr22_anno)
#' ftests <- ftest_orfs(psites %>% head(10000), chr22_anno, 
#'   n_cores=1)

ftest_orfs <- function(psites, anno, n_cores=1) {
  #
  orfs <- intersect(unique(psites$orf), names(anno$trspacecds))
  orfs <- anno$trspacecds[]
  psitecov <- psites %>%
    {
      x <- .
      out <- GenomicFeatures::mapToTranscripts(., orfs, ignore.strand = TRUE)
      out <- out[x$orf[out$xHits] == names(orfs)[out$transcriptsHits]]
      coverage(out)
    }
  #
  psitecov <- psitecov[sum(psitecov > 0) > 1]
  # now run multitaper tests on our data using multiple
  # cores if available
  message("running multitaper tests, this will be slow for a full dataset...")
  spec_tests <- psitecov %>% 
    parallel::mclapply(F = ftestvect, mc.cores = n_cores)
  # now format the output
  spec_test_df <- spec_tests %>%
    simplify2array() %>%
    t()
  spec_test_df <- spec_test_df %>% as.data.frame()
  spec_test_df$orf_id <- rownames(spec_test_df)
  rownames(spec_test_df) <- NULL
  colnames(spec_test_df) <- c("spec_coef", "p.value", "orf_id")
  spec_test_df$spec_coef <- spec_test_df$spec_coef %>% sqrt()
  orflens <- width(orfs[spec_test_df$orf_id])
  spec_test_df$spec_coef <- spec_test_df$spec_coef / orflens
  # put in NA values for things we couldn't test
  testdf <- tibble(orf_id = names(anno$cdsgrl)) %>% 
    left_join(spec_test_df, by = "orf_id")
  testdf
}



#' Run Multitaper tests on a set of ORFs
#'
#' @keywords Ribostan
#' @author Dermot Harnett, \email{dermot.p.harnett@gmail.com}
#'
#' @param psites GRanges object with psite information
#' @param anno  annotation object
#' @param remove  whether to remove non-peridic uORFs
#' @param n_cores  number of cores to use
#'
#' @details This function applies a multitaper test to
#' @return a numeric vector with the spectral coefficient at 0.333... and the pvalue for the test
#' @export
#' @examples
#' data(chr22_anno)
#' data(rpfs)
#' data(offsets_df)
#' psites <- get_psite_gr(rpfs, offsets_df, chr22_anno)
#' filteredanno <- periodicity_filter_uORFs(psites, chr22_anno)

periodicity_filter_uORFs <- function(psites, anno, remove=TRUE, n_cores=1){
  stopifnot(!is.null(anno$uORF))
  uORFs <- anno$uORF
  uORFs <- unique(names(uORFs[uORFs]))
  ftestdf <- ftest_orfs(psites, subset_annotation(anno, uORFs))
  #I guess we can just add to all elements
  if(is.null(mcols(anno$trspacecds)$spec_coef)) mcols(anno$trspacecds)$spec_coef <- NA
  if(is.null(mcols(anno$trspacecds)$p.value)) mcols(anno$trspacecds)$p.value <- NA
  mcols(anno$trspacecds[ftestdf$orf_id])<- ftestdf%>%select(-'orf_id')
  if(remove){
    periodic_uORFs <- ftestdf%>%filter(.data$p.value<0.05)%>%.$orf_id
    non_periodic_uORFs <- setdiff(uORFs, periodic_uORFs)
    orfs_to_keep <- setdiff(names(anno$trspacecds), non_periodic_uORFs)
    anno <- subset_annotation(anno, orfs_to_keep)
  }
  anno
}
