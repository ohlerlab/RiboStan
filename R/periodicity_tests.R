#' @importFrom multitaper dpss spec.mtm
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

ftestvect<-function(psit,k=24,bw=12){
  psit <- as.numeric(psit)
  slepians_values<-dpss(n=length(psit)%>%ifelse(.<25,50,.),k=k,nw=bw)
  # vals<-take_Fvals_spect(x = psit,n_tapers = k,time_bw = bw,slepians_values = sl)

  if(length(psit)<25){
      remain <- 50-length(psit)
      halfrmn <- as.integer(remain/2)
      psit<-c(rep(0,halfrmn),psit,rep(0,remain%%2+halfrmn))
  }
  #
  if(length(psit)<1024/2){padding<-1024}
  if(length(psit)>=1024/2){padding<-"default"}
  #
  resSpec1 <- spec.mtm(as.ts(psit), k=k, nw=bw, nFFT = padding, 
    centreWithSlepians = TRUE, Ftest = TRUE, 
    maxAdaptiveIterations = 100,returnZeroFreq=F,
    plot=F,dpssIN=slepians_values)
  psit
  #
  resSpec2<-dropFreqs(resSpec1,0.29,0.39)
  #
  closestfreqind <- which(abs((resSpec1$freq-(1/3)))==min(abs((resSpec1$freq-(1/3)))))
  # 
  freq_max_3nt<-resSpec1$freq[closestfreqind]
  Fmax_3nt<-resSpec1$mtm$Ftest[closestfreqind]
  spect_3nt<-resSpec1$spec[closestfreqind]
  return(c(Fmax_3nt,spect_3nt))

  pval <- pf(q=vals[1],df1=2,df2=(2*24)-2,lower.tail=F)
  return(c(vals[2],pval))
}


#' Run Multitaper Tests on All ORFs
#'
#' @keywords Ribostan
#' @author Dermot Harnett, \email{dermot.p.harnett@gmail.com}
#'
#' @param psites A GRanges object containing psites
#' @param anno annotation object
#'
#' @details This function applies a multitaper test to 
#' @return a numeric vector with the spectral coefficient at 0.333... and the pvalue for the test
#' @export

ftest_orfs<-function(psites, anno){
  #
  orfs <- anno$trspacecds[unique(psites$orf)]
  psitecov<-psites%>%
    {
        x <- .
        out <-GenomicFeatures::mapToTranscripts(.,orfs,ignore.strand=T)
        out <- out[x$orf[out$xHits]==names(orfs)[out$transcriptsHits]]
        coverage(out)
    }
  #
  ncores = detectCores()
  psitecov <- psitecov[sum(psitecov>0)>1]
  #now run multitaper tests on our data using multiple
  #cores if available
  message('running multitaper tests, this will be slow for a full dataset...')
  spec_tests <- psitecov%>%mclapply(ftestvect, mc.cores=ncores)
  #now format the output
  spec_test_df <- spec_tests%>%simplify2array%>%t
  spec_test_df <- spec_test_df%>%as.data.frame
  spec_test_df$orf_id <- rownames(spec_test_df)
  rownames(spec_test_df)<-NULL
  colnames(spec_test_df) <- c('spec_coef','p.value','orf_id')
  spec_test_df$spec_coef <- spec_test_df$spec_coef%>%sqrt
  orflens <- width(orfs[spec_test_df$orf_id])
  spec_test_df$spec_coef <- spec_test_df$spec_coef/orflens
  #put in NA values for things we couldn't test
  testdf <- tibble(orf_id = names(anno$cdsgrl))%>%left_join(spec_test_df)
}



