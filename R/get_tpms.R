library(tidyverse)
library(magrittr)
library(GenomicRanges)
library(Matrix)
id <- function(cov)match(cov,unique(cov))
#
ribobam <- '../cortexomics/pipeline/star_transcript/data/E13_ribo_1/E13_ribo_1.bam'
gtf <- rtracklayer::import('../cortexomics/ext_data/gencode.vM12.annotation.gtf')
tr2gid<-gtf%>%mcols%>%.[,c('gene_id','transcript_id')]%>%as.data.frame%>%distinct%>%filter(!is.na(transcript_id))%>%
	{setNames(.$gene_id,.$transcript_id)}
exonsgrl <- gtf%>%subset(type=='exon')%>%split(.,.$transcript_id)
#
ribogr <- GenomicAlignments::readGAlignments(ribobam,use.names=T)
seqlevels(ribogr)%<>%str_replace('\\|.*','')
names(ribogr) %<>% id
mcols(ribogr)$readlen <-  GenomicAlignments::qwidth(ribogr)
ribogr%<>%as("GenomicRanges")
#
if(seqnames(ribogr)%>%head(1)%>%str_detect('chr')){
	cov = ribogr%>%resize(1)%>%mapToTranscripts(exonsgrl[highcountcovtrs])
	cov$readlen = mcols(ribogr)$readlen[cov$xHits]
	cov$name = names(ribogr)[cov$xHits]
	cov%<>%subset(between(readlen,25,35))
	cov%<>%resize(1,'start')
} else{
	cov<-ribogr
	cov <- sort(cov)	
}
names(cov)%<>%id
library(Matrix)
spmat <- sparseMatrix(
	i=names(cov)%>%as.numeric,
	j=seqnames(cov)%>%id%>%as.vector,
	x=1
)
#now let's try the whole shebang in rstan
trlens <- exonsgrl%>%width%>%sum
sptrlens <- trlens[unique(cov@seqnames)%>%as.character]
fdata<-list(trlen=sptrlens)
fdata<-c(fdata,spmat%>%extract_sparse_parts)
fdata$trlen %<>% {./sum(.)}
fdata$TR <- spmat%>%ncol
fdata$R <- spmat%>%nrow
fdata$V <- fdata$w%>%length
fdata$Ulen <- fdata$u%>%length
fdata$classweights <- rep(1,fdata$R)
init = list(tpm=spmat%>%colSums%>%{./sum(.)})
#
eqtpm_mod <- rstan::stan_model(file='src/sparse_eq_tpm_opt.stan')
opt <- rstan::optimizing(eqtpm_mod,data=fdata,init=init,verbose=TRUE,iter=50)
#
tspmat <- t(spmat)
s_spmat <- summary(tspmat)
stopifnot(s_spmat$j%>%isSorted)
stopifnot(s_spmat$j%>%max%>%`>`(100e3))
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

#extract tpms from the stan object
tpms <- opt$par[opt$par%>%names%>%str_detect('tpm')]
#sample our data weighting by the tpms
matsample <- sample_cols_spmat(tspmat*tpms, return_mat=F)
#
iddf <- tibble(rind=1:length(cov),j=names(cov)%>%id,i=as.numeric(id(seqnames(cov))))
#
rinds <- iddf%>%inner_join(matsample,by=c('j','i'))%>%.$rind

sampcov <- cov[rinds]




#





#okay so I want this function to return a logical vector actually... 