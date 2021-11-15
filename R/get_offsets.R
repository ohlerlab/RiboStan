#see this script
#/fast/AG_Ohler/dharnet/cortexomics/src/Figures/Figure2_redo/1_p_alignment_redo.R:67
#for gr manip
#/fast/AG_Ohler/dharnet/cortexomics/src/Archive/Ribotransformer/exp_cds_4trans.R:94

################################################################################
########## Range manipulation
################################################################################
library(assertthat)
spl_mapFromTranscripts <- function(trspacegr,exons_grl){
  #
  exons_tr<-exons_grl%>%unlist%>%mapToTranscripts(exons_grl)%>%.[names(.)==seqnames(.)]
  ov <- findOverlaps(trspacegr,exons_tr)
  #
  trspacegr_spl <- suppressWarnings({trspacegr[queryHits(ov)]%>%pintersect(exons_tr[subjectHits(ov)])})
  genomic_trspacegr <- mapFromTranscripts(
  trspacegr_spl,
  # exons_tr[subjectHits(ov)]%>%split(.,seqnames(.))
  exons_grl
  )
  genomic_trspacegr$xHits <- queryHits(ov)[genomic_trspacegr$xHits]
  genomic_trspacegr
}
resize_grl_startfix<-function(grl,width){
  #what follows is some slightly black magic using S4 vectors
  #Integerlist which showings how much we'd need to trim that exon to get to to the desired transcript length
  trim =  cumsum(width(grl)) - width 
  #Where trim is greater than the exon width, we drop it
  drop = trim >=  width(grl)
  grl = grl[!drop]
  #vector showing location of the new 3' end of each transcript
  newends = cumsum(elementNROWS(grl))
  #vector with the amount we need to trim each new 3' end by
  endtrims=trim[IntegerList(as.list(elementNROWS(grl)))]@unlistData
  #finally, use these to trim
  grl@unlistData[newends] <- resize(grl@unlistData[newends], width(grl@unlistData[newends]) - endtrims  )
  grl
}
str_order_grl<-function(grl){BiocGenerics::order( start(grl)*(((strand(grl)!='-')+1)*2 -3) )}
sort_grl_st <- function(grl)grl[str_order_grl(grl),]
resize_grl_endfix <- function(grl,width){
  grl = invertStrand(grl)%>%sort_grl_st
  # 
  grl = resize_grl_startfix(grl,width)
  invertStrand(grl)%>%sort_grl_st
}

#' Perform a Ribo-seQC analysis
#'
#' This function loads annotation created by the prepare_annotation_files function, and analyzes a BAM file.
#' @param grl
#' @param width
#' @return resized grl
#' @details resize a grangeslist#' @import rmarkdown
#' @import rtracklayer
#' @import GenomicAlignments
#' @import BSgenome
#' @import GenomicFiles
#' @import devtools
#' @import reshape2
#' @import ggplot2
#' @import knitr
#' @import DT
#' @import gridExtra
#' @import ggpubr
#' @import viridis
#' @import Biostrings
#' @import GenomicFeatures
#' @import BiocGenerics
#' @import GenomicRanges
#' @export 

resize_grl <- function(grl,gwidth,fix='start',check=TRUE){
  stopifnot(all(gwidth>0))
  assert_that(all(all(diff(str_order_grl(grl))==1) ),msg = "grl needs to be 5'->3' sorted")
  if(fix=='start'){
    grl = resize_grl_startfix(grl,gwidth)
  }else if(fix=='end'){
    grl = resize_grl_endfix(grl,gwidth)
  }else if(fix=='center'){
    grlwidths = sum(width(grl)) 
    diffs = (gwidth - grlwidths)
    # 
    grl = resize_grl_startfix(grl,grlwidths + ceiling(diffs/2))
    grl = resize_grl_endfix(grl,grlwidths + diffs)
  }
  if(check){
    startstoolow <- any(start(grl)<=0)
    if(any(startstoolow)){
      stop(str_interp("${sum(startstoolow)} ranges extended below 1 .. e.g. ${head(which(startstoolow,1))}"))
    }
    grlseqs <- as.vector(unlist(use.names=F,seqnames(grl)[IntegerList(as.list(rep(1,length(grl))))]))
    endhighvect <- (GenomicRanges::end(grl)>GenomeInfoDb::seqlengths(grl)[grlseqs])
    endstoohigh <- any(endhighvect)
    if(any(endstoohigh)){
      stop(str_interp("${sum(endstoohigh)} ranges extended below above seqlength .. e.g. ${head(which(endstoohigh,1))}"))
    }
  }
  grl
}
trim_grl <- function(grl,bp,end='tp'){
  if(end=='tp'){
    resize_grl(grl,sum(width(grl)) - bp,fix='start')
  }else if(end=='fp'){
    resize_grl(grl,sum(width(grl)) - bp,fix='end')
  }else {
    stop("end should be fp or tp")
  }
}



################################################################################
########## Offsets
################################################################################
# trspacecds%<>%unlist
# trspacecds%>%width%>%`%%`(3)%>%table%>%{./sum(.)}

# trspacecds = trspacecds[is3bp]
# 
#

# 
# table((cds_prestop_st-cdsstarts)%%3)

get_orfcoords <- function(anno, fafile){

	library(Biostrings)
	library(Rsamtools)
	library(Rsamtools)
	#
	cdsgrl <- anno%>%subset(type=='CDS')%>%GenomicRanges::split(.,.$transcript_id)
	exonsgrl <- anno%>%subset(type=='exon')%>%GenomicRanges::split(.,.$transcript_id)
	#
	ribocovtrs <- names(cdsgrl)
	#
	trspacecds = GenomicFeatures::pmapToTranscripts(
		cdsgrl[ribocovtrs],
		exonsgrl[ribocovtrs])
	trspacecds <- unlist(trspacecds)
	stopifnot(all(names(trspacecds)==names(cdsgrl)))
	#chceck if the cds includes the stop codon
	fafile = Rsamtools::FaFile(fafile)
	cdsseqends = cdsgrl%>%
		sort_grl_st%>%
		resize_grl(3,'end')%>% 
		resize_grl(width(.)+3)%>%
		{.[width(.)%>%sum%>%`%%`(3)%>%`==`(0)];.}%>%
		sample(1000)%>%
		# resize(width(.)+3)%>%
		GenomicFeatures::extractTranscriptSeqs(x=fafile,.)%>%
		Biostrings::translate(.)
	cdsseqends = cdsseqends[Biostrings::nchar(cdsseqends)==2]
	#
	end_stop = table(subseq(cdsseqends,1,1))%>%sort%>%{./sum(.)}%>%tail(1)%>%`>`(0.5)
	end_plusone_stop = table(subseq(cdsseqends,2,2))%>%sort%>%{./sum(.)}%>%tail(1)%>%`>`(0.5)
	#set the cds_prestop_st as the 1st nuc of the last non-stop codon
	#
	cdsstarts = trspacecds%>%start%>%setNames(names(trspacecds))
	#
	if(end_plusone_stop){
		cds_prestop_st = trspacecds%>%end%>%`-`(2)%>%setNames(names(trspacecds))
	}else{
		cds_prestop_st = trspacecds%>%end%>%`-`(2-3)%>%setNames(names(trspacecds))
	}
	orf_coords <- data.frame(start=cdsstarts,prestop_st=cds_prestop_st)
	orf_coords
}


addphase <- function(gr,cdsstarts){
	gr$phase = unlist((start(gr) - cdsstarts[as.vector(seqnames(gr))])%%3)
	gr
}
string2onehot<-function(scol,acgt=c('A','C','T','G')){
	egval = rep(1,length(scol))
	acgt%>%
		vapply(function(lev){as.numeric(scol==lev)},egval)%>%
		matrix(ncol=4,dimnames=list(NULL,acgt))
}
#
dnaseq2onehot <- function(mat,pre){
	mat<-as.matrix(mat);
	nbp = ncol(mat)/2
	n=1
	posmats = lapply(1:ncol(mat),function(n){
		string2onehot(mat[,n,drop=FALSE])%>%
		set_colnames(paste0(pre,n-nbp-1,'_',colnames(.)))
	})%>%
	purrr::reduce(cbind)
	posmats
}
prop <- function(x,rnd=3) round(x/sum(x),rnd)
enddist <- function(gr){
	end(gr) - GenomeInfoDb::seqlengths(gr)[as.vector(seqnames(gr))]
}

width1grs <- function(gr){
	stopifnot(Negate(is.unsorted)(gr))
	isw1 <- width(gr)==1
	broad <- gr[!isw1]
	#vector of integers - 1,2,3 for range 1-3
	narrowstarts <- unlist(as(broad@ranges,'IntegerList'))
	narrow <- {GRanges(
			rep(seqnames(broad),width(broad)),
			IRanges(narrowstarts,w=1)
		)}
	mcols(narrow) <- mcols(broad)[rep(seq_along(broad),width(broad)),,drop=F]
	sort(c(gr[isw1],narrow))
}
# #TODO - make this count things lost gained rather than on off cds
# rln=29
# testtr=names(cdsstarts)[1]

#okay so this works if the 
get_offset_tbl <- function(ireads,readlens,cdsstarts, cds_prestop_st, ref_end='start'){
		readlens%<>%setNames(.,.)
		lapply(readlens,function(rln){
			#
			if(ref_end=='start'){
				tr_target_pos=cdsstarts_
			}else{
				tr_target_pos=cds_prestop_st_-2
			}
			ireads = ireads%>%
				subset(readlen==rln)%>%
				subset(seqnames%in%ribocovtrs)
			if(length(ireads)==0) return(NULL)
			#
			ireads = add_cor_offset(ireads,cdsstarts,cds_prestop_st)
			(cds_prestop_st-2-cdsstarts)%%3
			# return(length(ireads))
			#
			#
			allpos = expand.grid(phase=0:2,cor.offset=0:rln)%>%
				filter((cor.offset%%3)==0)%>%as.data.frame
			#
			stopifnot('score'%in%colnames(mcols(ireads)))
			frame_stat_df=ireads%>%
					subset(cor_offset<(rln-6))%>%
					subset(cor_offset>5)%>%
					{names(.)<-NULL;.}%>%
					as.data.frame%>%
					# mutate_at(vars(score),replace_na,0)%>%
					left_join(allpos,.)%>%
					group_by(phase,cor_offset)%>%
					tally(wt=score)%>%
					# tally()%>%
					# {(.$n*.$score)%>%sum}
					# group_by(phase)%>%
					ungroup%>%
					arrange(desc(cor_offset))%>%
					mutate(twind=n+lag(n)+lead(n))
					# group_by(cor_offset)%>%
					# filter(n()==3)%>%
					# mutate(twind=sum(n))
				#
			frame_stat_df%>%as.data.frame
			phasevect = 
				frame_stat_df%>%
				dplyr::slice(((which.max(twind)-1):(which.max(twind)+1)))%>%
				{setNames(.$cor_offset,.$phase)}
			message(capture.output(phasevect)%>%paste0(collapse='\n'))
			#
			# ireads%<>%{.$offset <- phasevect[as.character(.$phase)];.}
			# ireads%<>%{.$error=(.$cor_offset+.$phase)-(start(.)+.$offset);.}
			# ireads$cor_offset=NULL
			# list(ireads,frame_stat_df)
			frame_stat_df
	})
}
# #fake some data
# i=1
# cdsstarts[1]=101
# cds_prestop_st[1]=101+(3*30)
# # orfs <- 1:100%>%setNames(.,names(cdsstarts)[.])
# orfs <- 1%>%setNames(.,names(cdsstarts)[.])
# fakecov <- lapply(orfs,function(i){
# 	endbuff=6
# 	startbuff=6
# 	st = cdsstarts[[i]]
# 	nd = cds_prestop_st[[i]]+2
# 	#
# 	pos = seq(st-startbuff,nd+endbuff)
# 	#
# 	prob = (3 - (pos-st)%%3)/3
# 	#
# 	w = rbinom(length(prob),10,prob = prob)
# 	#
# 	w[startbuff+1] %<>% add(300) 
# 	w[1+startbuff+1] %<>% add(200) 
# 	w[2+startbuff+1] %<>% add(200) 
# 	w[3+startbuff+1] %<>% add(200) 
# 	w[4+startbuff+1] %<>% add(200) 
# 	w[length(w)-endbuff-2] %<>% add(300) 
# 	w[1+length(w)-endbuff-2] %<>% add(200) 
# 	w[2+length(w)-endbuff-2] %<>% add(200) 
# 	w[-1+length(w)-endbuff-2] %<>% add(200) 
# 	w[-2+length(w)-endbuff-2] %<>% add(200) 
# 	#
# 	ireads = data.frame(
# 		seqnames=names(cdsstarts)[[i]],
# 		start = pos,
# 		score=w
# 	)%>%mutate(end=start)
# 	#
# 	ireads
# })%>%bind_rows%>%mutate(readlen=29,end=start)%>%GRanges%>%
# 	GenomicRanges::shift(-12)%>%
# 	resize(29)
# ireads<-fakecov%>%resize(1)%>%coverage(weight='score')%>%GRanges%>%subset(score!=0)%>%width1grs
# ireads$readlen=29
#maybe we make this so it's aligning
add_cor_offset = function(ireads, cdsstarts, cds_prestop_st){
	ireads = ireads%>%
			as("GRanges")%>%
			addphase(cdsstarts)%>%
			{
				.$startoffset=cdsstarts[as.vector(seqnames(.))]+.$phase - start(.)
				.$endoffset=cds_prestop_st[as.vector(seqnames(.))]+.$phase - start(.)
				.$cor_offset = ifelse(
					(.$startoffset<(rln-5)) & (.$startoffset>5),
					.$startoffset,
					ifelse(
						(.$endoffset<(rln-5)) & (.$endoffset>5),
						.$endoffset,
						NA
					)
				)
				.
			}
	ireads	
}
# ireads <- ireads%>%add_cor_offset(cdsstarts, cds_prestop_st)
# #
# rl_offsets <- get_offset_tbl(ireads=ireads, 29, cdsstarts, cds_prestop_st)
# rl_offsets%<>%bind_rows(.id='length',.)
# #
# lphaseoffsetdf <- rl_offsets%>%
# 	group_by(length)%>%
# 	dplyr::slice(((which.max(twind)-1):(which.max(twind)+1)))

################################################################################
########## 
################################################################################

addphase <- function(gr,cdsstarts){
	gr$phase = unlist((start(gr) - cdsstarts[as.vector(seqnames(gr))])%%3)
	gr
}
add_cor_offset = function(ireads, cdsstarts, cds_prestop_st){
	ireads = ireads%>%
			as("GRanges")%>%
			addphase(cdsstarts)%>%
			{
				.$startoffset=cdsstarts[as.vector(seqnames(.))]+.$phase - start(.)
				.$endoffset=cds_prestop_st[as.vector(seqnames(.))]+.$phase - start(.)
				.$cor_offset = ifelse(
					(.$startoffset<(.$readlen-5)) & (.$startoffset>5),
					.$startoffset,
					ifelse(
						(.$endoffset<(.$readlen-5)) & (.$endoffset>5),
						.$endoffset,
						NA
					)
				)
				.
			}
	ireads	
}
# ireads<-sampled_cov_gr
get_incl_max_offsets<-function(sampled_cov_gr, anno, fafile){

	#
	cdsstarts <- anno$cdsstarts
	cds_prestop_st <- anno$cds_prestop_st
	#
	sampled_cov_gr <- sampled_cov_gr%>%
		add_cor_offset(cdsstarts, cds_prestop_st)
	#
	sampled_cov_gr$readlen%<>%as.numeric
	sampled_cov_gr$startoffset%>%table
	if(!'score'%in%colnames(mcols(sampled_cov_gr))) sampled_cov_gr$score<-1
	#data frame with numbers of reads gained at a specific offset
	gaindf = sampled_cov_gr%>%
		subset((startoffset<readlen-5) & (startoffset>5))%>%
		mcols(.)%>%.[,c('startoffset','score','readlen','phase')]%>%
		as.data.frame%>%
		group_by(phase, readlen, startoffset)%>%
		tally(wt=score)%>%
		arrange(startoffset)%>%
		group_by(phase, readlen)%>%
		mutate(gain=cumsum(n))%>%
		select(offset=startoffset,phase,gain)
	#data frame with numbers of reads lost at a specific offset
	lossdf = sampled_cov_gr%>%
		subset((.$endoffset<(readlen-3)) & (.$endoffset>3))%>%
		mcols(.)%>%.[,c('endoffset','score','readlen','phase')]%>%
		as.data.frame%>%
		group_by(phase, readlen, endoffset)%>%
		tally(wt=score)%>%
		arrange(endoffset)%>%
		group_by(phase, readlen)%>%
		mutate(loss=cumsum(n))%>%
		select(offset=endoffset,phase,loss)%>%
		mutate(offset=offset+3)
	#for a given offset, our gain is the start overlappers whose
	#offset is that or less, minus the stop guys whose offset is
	#less than that
	#
	netdf <- gaindf%>%
		full_join(lossdf, by=c('readlen', 'offset', 'phase'))%>%
		mutate(net=replace_na(gain,0)-replace_na(loss, 0))%>%
		group_by(readlen)%>%
		arrange(readlen,offset)%>%
		filter(!is.na(net))
	#
	readlens = sampled_cov_gr$readlen%>%unique
	#
	allpos = 
		map(readlens,function(rl){
			expand.grid(
				readlen=rl,
				phase=0:2,
				offset=2:(rl-3)
			)
		})%>%bind_rows%>%
		filter((offset%%3) ==0)
	netdf = allpos%>%left_join(netdf)
	netdf = netdf%>%
		mutate(net = replace_na(net,0))%>%
		mutate(twind = net+lag(net)+lead(net))
	netdf <- netdf%>%
		group_by(readlen)%>%
		# filter(readlen==17)%>%
		filter(max(net)>32)
	best_offsets <- netdf%>%
		filter(!is.na(twind))%>%
		group_by(readlen)%>%
		dplyr::slice((which.max(twind)-1):(which.max(twind)+1))
	#
	uniqueoffsetsfound = best_offsets%>%group_by(readlen, phase)%>%
		filter(n()>1)%>%nrow%>%`==`(0)
	stopifnot(uniqueoffsetsfound)
	#
	best_offsets
}


get_offsets <- function(sampled_cov_gr, anno, fafile, method='cds_incl'){
	# minreadlenfreq = sampled_cov_gr$readlen%>%table%>%{./sum(.)}%>%min
		# stopifnot(minreadlenfreq>0.05)
	if(method=='cds_incl'){
		best_offsets<-get_incl_max_offsets(sampled_cov_gr, anno, fafile)
	}else{
		stop('not yet implemented')
	}
	best_offsets
}