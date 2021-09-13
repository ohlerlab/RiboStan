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
#
ribogr <- GenomicAlignments::readGAlignments(ribobam,use.names=T)
seqlevels(ribogr)%<>%str_replace('\\|.*','')
names(ribogr) %<>% match(.,unique(.))
mcols(ribogr)$readlen <-  GenomicAlignments::qwidth(ribogr)
ribogr%<>%as("GenomicRanges")

seqnames(ribogr)[1:10000]%>%head(3)%>%id


tridrle <- seqnames(ribogr)%>%id
rdnmsrle <- Rle(names(ribogr))


inds=1:1e6
#
eqclasses=split(tridrle[inds],rdnmsrle[inds])
#
str_eqclasses<-eqclasses	%>%
	lapply(function(x){ out=x%>%	
		# {matrix(c(runLength(.),runValue(.)),byrow=T,nrow=2)%>%paste0(.,collapse=',')}
		{runValue(.)%>%paste0(.,collapse=',')}
		# stopifnot(is.character(out))
	})
#
eqclasscounts <- str_eqclasses%>%unlist%>%table

library(GenomicFeatures)


exonsgrl<-gtf%>%subset(type=='exon')%>%split(.,.$transcript_id)


ribogr%>%head(1000)%>%spl_mapFromTranscripts(exonsgrl)
ribogr%>%head(10e6)%>%GRanges%>%mapFromTranscripts(exonsgrl)


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
library(Matrix)
spmat <- sparseMatrix(
	i=names(cov)%>%as.numeric%>%id,
	j=seqnames(cov)%>%id%>%as.vector,
	x=1
)

spmat2 <- t(spmat)


for(i in 1:nrow(spmat2)){
	cat('.')
	spmat2 <- spmat2[,order(spmat2[i,])]
}


useqnames = unique(cov@seqnames)
tr2gid[useqnames]

#slowest
slowrle<-tr2gid[cov@seqnames%>%head(10e4)%>%as.vector]

#now, staying in Rle space
gid_rle <- Rle(
	tr2gid[seqnames(cov)%>%runValue%>%as.character],
	seqnames(cov)%>%runLength
)


library(Matrix)

g_spmat <- sparseMatrix(
	i=names(cov)%>%id,
	j=gid_rle%>%id%>%as.vector,
	x=1
)
gsharemat <- t(g_spmat) %*% g_spmat
library(igraph)
gsharemat[gsharemat!=0] <- 1
graphclust <- clusters(graph_from_graphnel(as(gsharemat,'graph')))

rs_spmat <- rowSums(spmat)
rstab <- table(rs_spmat)
rstab%>%cumsum%>%{./last(.)}

spmat[rs_spmat==1,]%>%{.%*%t(.)}





eq_list <- split(as.vector(id(seqnames(cov))),names(cov))
eq_list%>%head%>%table%>%as.array%>%str
summary(spmat)%>%group_by(i)%>%summarise(list(j))


trlist <- split(id(as.vector(seqnames(cov))),names(cov))%>%
	as('List')
trlist%>%sort
split(id(seqnames(cov)),names(cov))


scanbamparam <- ScanBamParam(what=c('qname','rname','pos','qwidth'))
bamscan<-scanBam(ribobams[1],param=scanbamparam)
bamscan <- bamscan[[1]]%>%lapply(Rle)

df <- DataFrame(r=bamscan[[1]],tr=bamscan[[2]])%>%
	subset(!is.na(tr))
df%>%
	
library(Matrix)


GenomicRanges()

library(Rcpp)
sourceCpp( "src/rowcount.cpp" )



rowCounts_3(Matrix::Matrix(sparse=T,1:24,ncol=3),1)
rowCounts_3(matrix(1:24,ncol=3),1)

tshare = t(spmat)%*%spmat
library(igraph)
tshare%>%dim
tshare <- 1*(tshare>0) 
graphclust <- clusters(graph_from_graphnel(as(tshare,'graph')))

cov



mmspmat <- spmat[rowSums(spmat)>1,]
mmspmat%*%t(mmspmat)

tspmat<-spmat%>%t
sum<-summary(spmat[rowSums(spmat)>1,])
eq <- sum%>%group_by(i)%>%summarise(paste(collapse=',',j))



rowCountsR <- function(x) {
  ## calculate hash
  m<-x
  h <- m %*% matrix(2^(0:(ncol(x)-1)), ncol=1)
  i <- 1:nrow(x)
  counts <- tabulate(h+1)
  counts[order(h[i])] <- counts
  list(counts=counts, idx=i)
}
rowCountsR(Matrix::Matrix(sparse=T,1:24,ncol=3))

segregate by rowSums I think

rs_rows <- spmat[rowSums(spmat)==2,]
rs_rows %*% t(rs_rows)



#


################################################################################
########## Test optimisation
################################################################################

#basic opt - works
tratio=1/3
fdata=list(
 map = 	rdat <- matrix(c(
	rep(c(1,0),250*tratio),
	rep(c(1,1),1000),
	rep(c(0,1),250)
	),byrow=T,ncol=2)
 ,
 trlen = c(1000+c(1000,1000))
)
fdata$trlen %<>% {./sum(.)}
fdata$TR <- fdata$map%>%ncol
fdata$R <- fdata$map%>%nrow
library(rstan)
tpm_mod <- rstan::stan_model(file='src/tpm_opt.stan')
opt <- rstan::optimizing(tpm_mod,data=fdata,algorithm="Newton",iter=10e4)
#this should get us a ratio of about 3
stopifnot(between(opt$par['tpm[1]']/opt$par['tpm[2]'],tratio-0.01,tratio+0.01))
opt


#now let's do this with equivalence classes
tratio=3
wfdata=list(
 map = 	rdat <- matrix(c(
	c(1,0),
	c(1,1),
	c(0,1)
	),byrow=T,ncol=2)
 ,
 classweights = matrix(c(tratio*250,1000,250)),
 trlen = c(1000+c(1000,1000))
)
wfdata$trlen %<>% {./sum(.)}
wfdata$TR <- wfdata$map%>%ncol
wfdata$R <- wfdata$map%>%nrow
library(rstan)
# for(i in 1:ncol(wfdata$map))wfdata$map[,i] = wfdata$map[,i]*wfdata$classweights
wmap = wfdata$map
for(i in 1:ncol(wmap))wmap[,i] = wfdata$map[,i]*wfdata$classweights
eqtpm_mod <- rstan::stan_model(file='src/eq_cl_tpm_opt.stan')
opt <- rstan::optimizing(eqtpm_mod,data=wfdata,algorithm="Newton")
opt

##WORKS!

#now with sparse matrix
tratio=3
wfdata=list(
 map = 	rdat <- matrix(c(
	c(1,0),
	c(1,1),
	c(0,1)
	),byrow=T,ncol=2)
 ,
 classweights = c(tratio*250,1000,250),
 trlen = c(1000+c(1000,1000))ls
)
wfdata<-c(wfdata,wfdata$map%>%extract_sparse_parts)
wfdata$trlen %<>% {./sum(.)}
wfdata$TR <- wfdata$map%>%ncol
wfdata$R <- wfdata$map%>%nrow
wfdata$V <- wfdata$w%>%length
wfdata$Ulen <- wfdata$u%>%length
library(rstan)
# for(i in 1:ncol(wfdata$map))wfdata$map[,i] = wfdata$map[,i]*wfdata$classweights
eqtpm_mod <- rstan::stan_model(file='src/sparse_eq_tpm_opt.stan')
opt <- rstan::optimizing(eqtpm_mod,data=wfdata,algorithm="Newton")
opt
#WORKS

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

opt <- rstan::optimizing(eqtpm_mod,data=fdata,init=init,verbose=TRUE,iter=50)

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
sample_cols_spmat <- function(spmat){
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
	passrand <- runif(spmat@x%>%length) < spmatnm@x 
	spmat@x <- xold
	spmat@x = spmat@x * passrand
	pass_summ = summary(spmat)
	pass_summ = pass_summ[pass_summ$x!=0,]
	pass_summ = as.data.frame(pass_summ)[diff(c(0,pass_summ$j))>0,]
	outmat = sparseMatrix(i=pass_summ$i,j=pass_summ$j,x=pass_summ$x,dims=dim(spmat))
	stopifnot(dim(outmat)==c(5,4))
	outmat	
}
#spmat


