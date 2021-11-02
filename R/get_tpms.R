library(tidyverse)
library(magrittr)
library(GenomicRanges)
library(Matrix)

id <- function(cov) BiocGenerics::match(cov,unique(cov))



process_ribogr <- function(ribobam, strip_seqnames=TRUE){
	#
	require(Rsamtools)
	require(GenomicAlignments)
	bparam <- ScanBamParam(simpleCigar = TRUE,scanBamFlag(isUnmappedQuery = FALSE, 
		isMinusStrand = FALSE, ))
	ribogr <- GenomicAlignments::readGAlignments(ribobam,use.names=T,param=bparam)
	#
	wfilt = qwidth(ribogr)<=35 & (20<=qwidth(ribogr))
	ribogr = ribogr[wfilt]
	# this is the number get_bamdf gets in py 4628644
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
	# names(cov)%<>%id
	cov	
}

get_cds_reads<-function(cov, anno){
	trspacecds <- anno$trspacecds
	sharedseqnames <- unique(seqnames(trspacecds))%>%intersect(seqnames(cov))
	cov <-  cov%>%keepSeqlevels(sharedseqnames,pruning='coarse')
	seqlevels(cov) <- seqlevels(trspacecds)
	seqinfo(cov) <- seqinfo(trspacecds)
	cov <- cov %>% subsetByOverlaps(trspacecds)
}
get_readgr<-function(ribobam, anno, strip_seqnames=TRUE){
	bamseqnames <- seqinfo(Rsamtools::BamFile(ribobam))@seqnames
	if(strip_seqnames) bamseqnames%<>%str_replace('\\|.*','')
	stopifnot(mean(unique(seqnames(anno$trspacecds)) %in% bamseqnames)>.9)
	cov <- process_ribogr(ribobam)
	cov <- get_cds_reads(cov, anno)
	#4542346 is the number cds filtering gets us down to in py
	names(cov)%<>%id
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

	spmat <- Matrix::sparseMatrix(
		i=names(cov)%>%as.numeric,
		j=seqnames(cov)%>%id%>%as.vector,
		x=1
	)
	colnames(spmat) <- seqnames(cov)%>%unique
	spmat = spmat / Matrix::rowSums(spmat)
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
	trlens <- anno$trspacecds%>%width%>%
		setNames(names(anno$trspacecds))
	#now let's try the whole shebang in rstan
	sptrlens <- trlens[colnames(spmat)]
	fdata<-list(trlen=sptrlens)
	fdata<-c(fdata,spmat%>%rstan::extract_sparse_parts(.))
	fdata$trlen %<>% {./sum(.)}
	fdata$TR <- spmat%>%ncol
	fdata$R <- spmat%>%nrow
	fdata$V <- fdata$w%>%length
	fdata$Ulen <- fdata$u%>%length
	fdata$classweights <- rep(1,fdata$R)
	init = list(tpm=spmat%>%{Matrix::colSums(.)}%>%`/`(fdata$trlen)%>%{./sum(.)})
	#
	modelcode = '
		data {
		  int TR;// number of TRs
		  int R;// number of reads
		  int V;
		  int Ulen;
		  vector [V] w ;
		  int  v [V];
		  int   u [Ulen];
		  vector [R] classweights;
		}
		parameters {
		  simplex [TR] n;
		}
		model {
		    target += log(
		      csr_matrix_times_vector(R, TR, w, v, u, n)
		    ).* classweights;
		}
	'
	eqtpm_mod <- rstan::stan_model(model_code = modelcode)
	opt <- rstan::optimizing(
		eqtpm_mod,
		data=fdata,
		init=init,
		verbose=FALSE,
		iter=iternum)
	opt$seqnames <- colnames(spmat)
	opt$trlen <- fdata$trlen
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
	if(tpm_opt$par%>%names%>%str_detect('^tpm\\[\\d+\\]$')%>%any){
		tpmpars <- tpm_opt$par%>%names%>%str_subset('^tpm\\[\\d+\\]$')
		tpms <- tpm_opt$par[tpmpars]
		names(tpms) <- tpm_opt$seqnames
		tpms <- tpms * 1e6
		tpms	
	}else{
		tpmpars <- tpm_opt$par%>%names%>%str_subset('^n\\[\\d+\\]$')
		tpms <- tpm_opt$par[tpmpars]
		names(tpms) <- tpm_opt$seqnames
		tpms = tpms / tpm_opt$trlen
		tpms = tpms / sum(tpms)
		tpms <- tpms * 1e6
		tpms	
	}
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


gene_level_expr <- function(tpms, anno){
	trgiddf = anno$trgiddf
		#
	gn_expr=left_join(
		trgiddf,
		enframe(tpms, 'transcript_id', 'tpm'),
		by='transcript_id')
	#	
	gn_expr%<>%group_by(gene_id)%>%summarise(expr = sum(replace_na(tpm,0)))
	gn_expr%>%select(gene_id, expr)
}

get_ribofasta_anno <- function(ribofasta){
	Rsamtools::indexFa(ribofasta)
	faheadernames <- seqinfo(Rsamtools::FaFile(ribofasta))
	faheaddf=seqnames(faheadernames)%>%as.vector%>%str_split_fixed('\\|',10)
	anno = list()
	anno$trspacecds = GRanges(faheaddf[,1],faheaddf[,9]%>%str_extract('\\d+\\-\\d+')%>%
		str_split_fixed('-',2)%>%set_colnames(c('start','end'))%>%
		apply(2,as.numeric)%>%as.data.frame%>%{IRanges(start=.$start,end=.$end)})%>%
		setNames(.,as.character(seqnames(.)))
	anno$trgiddf <- tibble(transcript_id = faheaddf[,1], gene_id = faheaddf[,2])
	anno = c(anno,
		list(cdsstarts = anno$trspacecds%>%start%>%
			setNames(names(anno$trspacecds)),
		cds_prestop_st = anno$trspacecds%>%end%>%`-`(2)%>%
			setNames(names(anno$trspacecds))	
	))
	anno
}

get_exprfile <- function(ribobam, ribofasta, outfile){
	#	
	anno <- get_ribofasta_anno(ribofasta)
	#
	cov <- get_readgr(ribobam, anno)
	#
	spmat <- get_read_spmat(cov)
	#
	tpm_opt <- optimize_tpms(spmat, anno, iternum=100)
	#
	tpms <- get_tpms(tpm_opt)
	#
	n_reads = nrow(spmat)
	counts = n_reads * (tpms / 1e6)
	counts = enframe(counts, 'Name', 'NumReads')
	tpmdf = enframe(tpms, 'Name', 'TPM')
	#
	cdslens <- anno$trspacecds%>%width%>%setNames(names(anno$trspacecds))%>%
		enframe('Name','Length')%>%mutate(EffectiveLength=Length)
	output <- cdslens%>%left_join(tpmdf)%>%left_join(counts)
	#
	output <- output %>% select(Name, Length, EffectiveLength, TPM, NumReads)
	output%>%write_tsv(outfile)
}


