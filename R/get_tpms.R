# RiboseQC, a comprehensive Ribo-seq quality control tool
#
# Authors:
# Lorenzo Calviello (calviello.l.bio@gmail.com)
# Dominique Sydow (dominique.sydow@posteo.de)
# Dermott Harnett (Dermot.Harnett@mdc-berlin.de)
# Uwe Ohler (Uwe.Ohler@mdc-berlin.de)
#
# This software is free software: you can redistribute it and/or
# modify it under the terms of the GNU General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.
#
# This software is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this software. If not, see
# <http://www.gnu.org/licenses/>.

#' @import GenomicFeatures
#' @import GenomicFiles
NULL


library(tidyverse)
library(magrittr)
library(GenomicRanges)
library(Matrix)

id <- function(cov) BiocGenerics::match(cov,unique(cov))


#' Select those read lengths which make up 95% of the reads
#'
#' This function takes in a vector of the lengths of each read, and determines which read lengths
#' should be included in order to filter out the 2.5% shortest/longest reads.
#'
#' @keywords Ribostan
#' @author Dermot Harnett, \email{dermot.p.harnett@gmail.com}
#'
#' @param readlens Numeric; numeric vector
#' @return a numeric integer vector, such that selecting reads with these lengths will preserve at least 95%
#' of reads

get_readlens <- function(readlens){
	readlens%>%as.numeric%>%
	table%>%cumsum%>%{./max(.)}%>%
	keep(~.>0.025)%>%
	{keep(.,~1-.>0.025)}%>%
	names%>%
	as.numeric
}


#' Read a bam file containing Ribosomal Footprints
#'
#' This function 
#'
#' @keywords Ribostan
#' @author Dermot Harnett, \email{dermot.p.harnett@gmail.com}
#'
#' @param ribobam String; full path to html report file.
#' @param strip_seqnames  whether the function should remove all text after the first '|', useful if aligning to gencode fastas#
#' Defaults to \code{TRUE}
#' @return This function returns a granges object with the read names in the names slot and a metadata column
#' denoting readlength.
#'
#' @details This function reads in a bam file containing ribosomal footprints. It uses only reads without splice sites
#' and reads which align to the positive strand, as it's designed to work on transcriptomic alignments.
#' 
#' @seealso \code{\link{get_cds_reads}}, \code{\link{get_readlens}}

read_ribobam <- function(ribobam, strip_seqnames=TRUE){
	#
	require(Rsamtools)
	require(GenomicAlignments)
	bparam <- ScanBamParam(simpleCigar = TRUE,scanBamFlag(isUnmappedQuery = FALSE, 
		isMinusStrand = FALSE, ))
	ribogr <- GenomicAlignments::readGAlignments(ribobam,use.names=T,param=bparam)
	#
	readlens <- get_readlens(qwidth(ribobam))
	wfilt = qwidth(ribogr)<=max(readlens) & (min(readlens)<=qwidth(ribogr))
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
		cov = ribogr%>%resize(1)%>%mapToTranscripts(exonsgrl)
		cov$readlen = mcols(ribogr)$readlen[cov$xHits]
		cov$name = names(ribogr)[cov$xHits]
		cov%<>%subset(between(readlen, min(readlens), max(readlens)))
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

#' Filter an RPF GR for overlap with coding sequences
#'
#' @keywords Ribostan
#' @author Dermot Harnett, \email{dermot.p.harnett@gmail.com}
#'
#' @param cov GRanges; A granges object with 
#' @param strip_seqnames  whether the function should remove all text after the first '|', useful if aligning to gencode fastas#
#' Defaults to \code{TRUE}
#' @return This function returns a granges object with the read names in the names slot and a metadata column
#' denoting readlength.
#'
#' @details This function reads in a bam file containing ribosomal footprints. It uses only reads without splice sites
#' and reads which align to the positive strand, as it's designed to work on transcriptomic alignments.
#' 
#' @seealso \code{\link{get_cds_reads}}, \code{\link{get_readlens}}

get_cds_reads<-function(cov, anno){
	trspacecds <- anno$trspacecds
	sharedseqnames <- unique(seqnames(trspacecds))%>%
		intersect(seqnames(cov))%>%
		unlist
	cov <-  cov%>%keepSeqlevels(sharedseqnames,pruning='coarse')
	seqlevels(cov) <- seqlevels(trspacecds)
	seqinfo(cov) <- seqinfo(trspacecds)
	cov <- cov %>% subsetByOverlaps(trspacecds)
}


#' Read and filter a bam file of RPFs
#'
#' @keywords Ribostan
#' @author Dermot Harnett, \email{dermot.p.harnett@gmail.com}
#'
#' @param ribobam GRanges; A bam file with RPFs
#' @param ribobam GRanges; A coverage
#' @param strip_seqnames  whether the function should remove all text after the first '|', useful if aligning to gencode fastas#
#' Defaults to \code{TRUE}
#' @return This function returns a granges object with the read names in the names slot and a metadata column
#' denoting readlength.
#'
#' @details This function reads in a bam file containing ribosomal footprints. It uses only reads without splice sites
#' and reads which align to the positive strand, as it's designed to work on transcriptomic alignments.
#' 
#' @seealso \code{\link{get_cds_reads}}, \code{\link{get_readlens}}
#' @details This function creates the html report visualizing the RiboseQC analysis data. \cr \cr
#' Input are two lists of the same length: \cr \cr
#' a) \code{input_files}: list of full paths to one or multiple input files
#' (RiboseQC analysis files generated with \code{RiboseQC_analysis}) and \cr \cr
#' b) \code{input_sample_names}: list of corresponding names describing the file content in max. 5 characters
#' (these are used as names in the report). \cr \cr
#' For the report, a RMarkdown file is rendered as html document, saved as \code{output_file}. \cr \cr
#' Additionally, all figures in the report are saved as PDF figures
#' in an extra folder in the same directory as the report html file. \cr \cr
#' Example: \cr
#' \code{output_file <- "\\mydir\\myreport.html"} will generate
#' the html report \code{\\mydir\\myreport.html} and
#' the folder \code{\\mydir\\myreport_plots\\} for the RDS object files to be stored in.
#'
#' @seealso code{\\mydir\\myreport_plots\\}
#' 
#' @examples
#' \dontrun{
#' 
#' }
#' @export

get_readgr<-function(ribobam, anno, strip_seqnames=TRUE){
	bamseqnames <- seqinfo(Rsamtools::BamFile(ribobam))@seqnames
	if(strip_seqnames) bamseqnames%<>%str_replace('\\|.*','')
	stopifnot(mean(unlist(unique(seqnames(anno$trspacecds))) %in% bamseqnames)>.9)
	cov <- read_ribobam(ribobam)
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
		unlist%>%
		setNames(names(anno$trspacecds))
	#now let's try the whole shebang in rstan
	sptrlens <- trlens[colnames(spmat)]
	fdata<-list(nonorm_trlen=sptrlens)
	fdata<-c(fdata,spmat%>%rstan::extract_sparse_parts(.))
	fdata$trlen <- fdata$nonorm_trlen%>%{./sum(.)}
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
	opt$trlen <- fdata$nonorm_trlen[opt$seqnames]
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


