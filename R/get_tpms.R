library(tidyverse)
library(magrittr)
library(GenomicRanges)
library(Matrix)

################################################################################
########## 
################################################################################
filter_anno <- function(anno){

	filt_anno <- anno

	{
		require(Biostrings)
		#find which cds are multiples of 3bp
		cdsgrl <- filt_anno%>%subset(type=='CDS')%>%split(.,.$transcript_id)
		is3bp = cdsgrl%>%width%>%sum%>%`%%`(3)%>%`==`(0)	
		ribocovtrs = names(cdsgrl)[is3bp]
		#subset cds and anno with these
		cdsgrl <- cdsgrl[ribocovtrs]
		filt_anno = filt_anno%>%
			subset(type!='CDS')%>%
			subset(transcript_id%in%ribocovtrs)
		#chceck if the cds includes the stop codon
		fafileob = Rsamtools::FaFile(fafile)
		cdsgrl%<>%sort_grl_st
		cdsseqends = cdsgrl%>%
			resize_grl(3,'end')%>% 
			resize_grl(sum(width(.))+3)%>%
			GenomicFeatures::extractTranscriptSeqs(x=fafileob,.)
		#some sequences have spaces (ends of chrs i think)	
		filterchars = cdsseqends%>%str_detect('[^ATCG]')
		cdsseqends[filterchars] = 'AAAAAA'
		cdsseqends %<>% translate
		stopifnot(cdsseqends%>%nchar%>%is_in(2))
		# cdsseqends = cdsseqends[nchar(cdsseqends)==2]
		#
		end_stop = table(subseq(cdsseqends,1,1))%>%sort%>%
			{./sum(.)}%>%.['*']%>%`>`(0.5)
		end_plusone_stop = table(subseq(cdsseqends,2,2))%>%sort%>%
			{./sum(.)}%>%.['*']%>%`>`(0.5)
		stopifnot(end_stop|end_plusone_stop)
		#
		endseq = if(end_plusone_stop){ subseq(cdsseqends,2,2) }else{
		 subseq(cdsseqends,1,1)
		}
		#set the cds_prestop_st as the 1st nuc of the last non-stop codon
		#
		exonsgrl <- filt_anno%>%
			subset(type=='exon')%>%
			split(.,.$transcript_id)%>%
			.[ribocovtrs]
		stopifnot(all(ribocovtrs==names(cdsgrl)))
		stopifnot(all(ribocovtrs==names(exonsgrl)))
		trspacecds = GenomicFeatures::pmapToTranscripts(
			cdsgrl[ribocovtrs],
			exonsgrl[ribocovtrs])
		stopifnot(all(ribocovtrs==names(trspacecds)))		
		trspacecds <- trspacecds%>%unlist
		stopifnot(all(ribocovtrs==names(trspacecds)))		


		#
		cdsstarts = trspacecds%>%start%>%setNames(names(trspacecds))
		#
		if(end_stop) cdsgrl <- cdsgrl%>%resize_grl(sum(width(.))-3,'start')
		stopifnot(all(ribocovtrs==names(cdsgrl)))
		filt_anno = c(filt_anno,unlist(cdsgrl))
		#now only those which have M at the start and '*' at the end
		cdsseqstarts = cdsgrl%>%
			sort_grl_st%>%
			resize_grl(3,'start')%>% 
			GenomicFeatures::extractTranscriptSeqs(x=fafileob,.)%>%
			Biostrings::translate(.)
		#	
		#	
		hasMstart = cdsseqstarts[ribocovtrs]=='M'
		hasStop = endseq[ribocovtrs]=='*'

		ribocovtrs = ribocovtrs[hasMstart & hasStop]
		outanno = list(
			ribocovtrs=ribocovtrs,
			trspacecds=trspacecds[ribocovtrs])
		outanno = c(outanno,
			list(cdsstarts = outanno$trspacecds%>%start%>%
				setNames(names(outanno$trspacecds)),
			cds_prestop_st = outanno$trspacecds%>%end%>%`-`(2)%>%
				setNames(names(outanno$trspacecds)),
			anno = filt_anno%>%subset(transcript_id%in%ribocovtrs)		
		))
	}
	return(outanno)
}


id <- function(cov)match(cov,unique(cov))



process_ribogr <- function(ribobam, strip_seqnames=TRUE){
	#
	ribogr <- GenomicAlignments::readGAlignments(ribobam,use.names=T)
	#
	names(ribogr) %<>% id
	mcols(ribogr)$readlen <-  GenomicAlignments::qwidth(ribogr)
	ribogr%<>%as("GenomicRanges")
	#strip seqnames for when bam files contains a large fasta header
	if(strip_seqnames){
		seqlevels(ribogr)%<>%str_replace('\\|.*','')
	}
	#load genomic or transcriptomic bam
	if(seqnames(ribogr)%>%head(1)%>%str_detect('chr')){
		cov = ribogr%>%resize(1)%>%mapToTranscripts(exonsgrl[highcountcovtrs])
		cov$readlen = mcols(ribogr)$readlen[cov$xHits]
		cov$name = names(ribogr)[cov$xHits]
		cov%<>%subset(between(readlen,25,35))
		cov%<>%resize(1,'start')
		cov <- sort(cov)
	} else{
		cov<-ribogr
		cov <- sort(cov)	
	}
	#name the reads with integers as they appear in the sorted object
	names(cov)%<>%id
	cov	
}


get_cds_reads<-function(cov, anno){
	cdsgrl <- anno$anno%>%subset(type=='CDS')%>%split(.,.$transcript_id)
	exonsgrl <- anno$anno%>%subset(type=='exon')%>%split(.,.$transcript_id)
	#
	ribocovtrs <- names(cdsgrl)
	#
	trspacecds <- GenomicFeatures::pmapToTranscripts(
		cdsgrl[ribocovtrs],exonsgrl[ribocovtrs])
	trspacecds <- trspacecds %>% unlist
	valid_cds <-  names(trspacecds)==names(cdsgrl)
	stopifnot(valid_cds)
	cov <- cov %>% subsetByOverlaps(trspacecds)
	names(cov)%<>%id
	cov
}
	

get_readgr<-function(ribobam, anno, strip_seqnames=TRUE){
	bamseqnames <- seqinfo(Rsamtools::BamFile(ribobam))@seqnames
	if(strip_seqnames) bamseqnames%<>%str_replace('\\|.*','')
	stopifnot(mean(unique(anno$anno$transcript_id) %in% bamseqnames)>.9)
	cov <- process_ribogr(ribobam)
	cov <- get_cds_reads(cov, anno)
	cov
}


get_cds_profile <- function(cov,trspacecds){
	cdsstarts <- setNames(start(trspacecds),names(trspacecds))
	cov%>%subset(readlen=='29')
	cov$cdsstart <- cdsstarts[as.character(seqnames(cov))]
	cov$cdspos <- start(cov)-cov$cdsstart

	cdsends <- setNames(end(trspacecds),names(trspacecds))
	cov$cdsend <- cdsends[as.character(seqnames(cov))]
	cov$endcdspos <- start(cov)-cov$cdsend

	cdspos <- cov%>%
		subset(cdspos%>%between(-300,300))%>%
		subset(endcdspos%>%`<`(-300))%>%
		.$cdspos

		library(txtplot)

	cdspos%>%table%>%.[as.character(-30:30)]%>%replace_na(0)%>%
		txtplot

}


#' get_read_spmat
#'
#' get optimal tpms using stan
#'
#' @param cov GRanges object of 
#' @return A matrix of the infile
#' @export

get_read_spmat <- function(cov){
	require(Matrix)

	spmat <- sparseMatrix(
		i=names(cov)%>%as.numeric,
		j=seqnames(cov)%>%id%>%as.vector,
		x=1
	)
	colnames(spmat) <- seqnames(cov)%>%unique
	spmat
}



#' optimize_tpms
#'
#' get optimal tpms using stan
#'
#' @param spmat a sparse numeric matrix
#' @param trlens a sparse numeric matrix
#' @return A matrix of the infile
#' @export
optimize_tpms<-function(spmat, anno, iternum=500){
	require(rstan)
	#
	trlens <- anno%>%subset(type=='CDS')%>%
		split(.,.$transcript_id)%>%
		.[colnames(spmat)]%>%
		width%>%sum
	trlens <- anno%>%subset(type=='exon')%>%
		split(.,.$transcript_id)%>%
		.[colnames(spmat)]%>%
		width%>%sum
	#now let's try the whole shebang in rstan
	sptrlens <- trlens[unique(cov@seqnames)%>%as.character]
	fdata<-list(trlen=sptrlens)
	fdata<-c(fdata,spmat%>%rstan::extract_sparse_parts(.))
	fdata$trlen %<>% {./sum(.)}
	fdata$TR <- spmat%>%ncol
	fdata$R <- spmat%>%nrow
	fdata$V <- fdata$w%>%length
	fdata$Ulen <- fdata$u%>%length
	fdata$classweights <- rep(1,fdata$R)
	init = list(tpm=spmat%>%colSums%>%{./sum(.)})
	spmat%>%colSums%>%max
	#
	eqtpm_mod <- rstan::stan_model(file='src/sparse_eq_tpm_opt.stan')
	opt <- rstan::optimizing(
		eqtpm_mod,
		data=fdata,
		# init=init,
		verbose=TRUE,
		iter=iternum)
	opt$seqnames <- colnames(spmat)
	opt
}


#
# tspmat <- t(spmat)
# s_spmat <- summary(tspmat)
# stopifnot(s_spmat$j%>%isSorted)
# stopifnot(s_spmat$j%>%max%>%`>`(100e3))

#TODO now write the code for sampling the matrix in a weighted way once
#we have tpms. 
# spmat<-Matrix(c(0,1,0,0 
# 								,1,0,0,1, 
# 								 0,3,0,1,
# 								 0,1,1,1,
# 								 0,1,1,50),nrow=5,ncol=4,byrow=T,sparse=T)
# spmat@i#0 indexed row
# spmat@p#initial element for each column, 0 based.
# spmat@x
#take in a sparse matrix - sample one value from each of it's columns,
#weighted by the value of those columns

#' sample_cols_spmat
#'
#' This function takes in a sparse matrix, and then samples a single value
#' from each column, with sampling weights within each column equal to
#' the columns values
#'
#' @param spmat a sparse numeric matrix
#' @param infile 
#' @return A matrix of the infile
#' @export


sample_cols_spmat <- function(spmat,return_mat=T){
	#	
	#so I think p gives us element which is the final one for each 
	matsum <- summary(spmat)
	cs <- cumsum(spmat@x)
	colsums <- colSums(spmat)
	pts <- spmat@p
	prevcs <- c(0,cs[pts%>%tail(-1)])[matsum$j]
	cs = cs - prevcs
	xold <- spmat@x
	spmat@x <- cs
	spmat
	spmatnm = t(t(spmat) / colsums)
	nvals = spmat@x%>%length
	passrand <- runif(nvals) < spmatnm@x 
	# firstpass = diff(c(0,passrand))==1
	# logmat = spmat>0
	# logmat@x = passrand
	# (spmat*logmat)%>%summary%>%subset(x!=0)%>%nrow
	#
	spmat@x <- xold
	spmat@x = spmat@x * passrand
	pass_summ = summary(spmat)%>%as.data.frame
	# pass_summ$ind = 1:nrow(pass_summ)
	pass_summ = pass_summ[pass_summ$x!=0,]
	pass_summ = pass_summ[diff(c(0,pass_summ$j))>0,]
	#returnmatrix or the summary
	if(return_mat){
		outmat = sparseMatrix(i=pass_summ$i,j=pass_summ$j,x=pass_summ$x,dims=dim(spmat))
		stopifnot(dim(outmat)==dim(spmat))		
		outmat
	}else{
		pass_summ
	}
}

################################################################################

#' get_tpms
#'
#' This function takes in a sparse matrix, and then samples a single value
#' from each column, with sampling weights within each column equal to
#' the columns values
#'
#' @param tpm_opt a sparse numeric matrix
#' @return a vector of normalized footprint densities
#' @export


get_tpms <- function(tpm_opt){
	tpmpars <- tpm_opt$par%>%names%>%str_subset('^tpm\\[\\d+\\]$')
	tpms <- tpm_opt$par[tpmpars]
	names(tpms) <- tpm_opt$seqnames
	tpms <- tpms * 1e6
	tpms
}

#' sample_cov_gr
#'
#' This function takes in a coverage GR and takes one out of each multimap
#' weighting according to the tpms 
#' 
#' @param tpm_opt a sparse numeric matrix
#' @return a vector of normalized footprint densities
#' @export


sample_cov_gr <- function(cov, tpms){
	spmat <- get_read_spmat(cov)
	#
	matsample <- sample_cols_spmat(t(spmat)*tpms, return_mat=F)
	#
	iddf <- tibble(rind=1:length(cov),j=names(cov)%>%id,i=as.numeric(id(seqnames(cov))))
	#
	rinds <- iddf%>%inner_join(matsample,by=c('j','i'))%>%.$rind
	#
	sampcov <- cov[rinds]
	#
	sampcov
}

