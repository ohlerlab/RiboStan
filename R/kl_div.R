
# if(!file.exists(here('data/cdsfpcovlist.rds'))){
# 	cdsfpcovlist <- ribobams%>%mclapply(mc.cores=8,function(ribobam){
# 		ribogr <- GenomicAlignments::readGAlignments(ribobam)
# 		mcols(ribogr)$readlen <-  GenomicAlignments::qwidth(ribogr)
# 		ribogr%<>%as("GenomicRanges")
# 		cov = ribogr%>%resize(1)%>%mapToTranscripts(exonsgrl[names(innercdspos)])
# 		cov$readlen = mcols(ribogr)$readlen[cov$xHits]
# 		cov%<>%subset(between(readlen,LOWREADLIM,HIGHREADLIM))
# 		cov%<>%resize(1,'start')
# 		cdscov <- cov%>%mapToTranscripts(innercdspos)
# 		cdscov$readlen <- cov$readlen[cdscov$xHits]
# 		cdscov$xHits <- NULL
# 		cdscov$transcriptsHits<-NULL
# 		split(cdscov,cdscov$readlen)%>%
# 			lapply(coverage)
# 	})
# 	cdsfpcovlist%<>%setNames(samples)
# 	saveRDS(cdsfpcovlist,here('data/cdsfpcovlist.rds'))
# }else{
# 	cdsfpcovlist<-readRDS(here('data/cdsfpcovlist.rds'))
# }



get_codon_matches <- function(codon, anno,fafileob, n_start_buff,n_end_buff,n_wind_l_ext,n_wind_r_ext){
		exonsgrl=anno$exonsgrl
		trspacecds=anno$trspacecds
		#define inner cds - cds but with start and end trimmed off.
		innercds = trspacecds%>%subset(width>(3+n_start_buff+n_end_buff))%>%
					resize(width(.)-n_start_buff,'end')%>%
					resize(width(.)-n_end_buff,'start')
		exonseq = exonsgrl%>%extractTranscriptSeqs(x=fafileob)
		#get matches to that codon (any frame)
		codmatches<-vmatchPattern(pattern=codon,exonseq)#exclude the start ccodon
		matchgr<-codmatches%>%unlist%>%GRanges(names(.),.)
		#select only the in frame codons
		matchgr$cdspos = start(matchgr) - start(trspacecds[as.vector(seqnames(matchgr))])
		matchgr%<>%subset(cdspos %%3 == 0)
		#put seqlengths on this object
		seqlengths(matchgr) = exonsgrl%>%width%>%.[seqlevels(matchgr)]%>%sum	
		#use only matches in the inner cds 
		matchgr = matchgr%>%subsetByOverlaps(innercds)
		#expand the windows around these codons
		codmatchwindows<-matchgr%>%
			resize(n_wind_l_ext,'end')%>%
			resize(n_wind_r_ext,'start')
		#require these to be inside our cds
		codmatchwindows <- codmatchwindows[!is_out_of_bounds(codmatchwindows)]
		codmatchwindows%<>%subsetByOverlaps(innercds)
		innercdspos <- innercds%>%{strand(.)<-'+';.}
		#lift them to inner cds coordinates
		icds_codons <- allcodlist%>%mapToTranscripts(innercdspos)
		mcols(icds_codons)=NULL
		icds_codons
	}

get_cds_codons <- function(anno, fafileob, 
	n_wind_l_ext = 45,
	n_wind_r_ext = 9,
	n_start_buff=60,
	n_end_buff=60,
	allcodons=getGeneticCode())
{
	allcodons <- allcodons%>%setNames(names(allcodons))
	cds_codons <- lapply(
		allcodons,
		anno = anno,
		fafileob = fafileob,
		n_start_buff=n_start_buff, n_end_buff=n_end_buff,
		n_wind_l_ext=n_wind_l_ext, n_wind_l_ext=n_wind_l_ext,
		F=get_codon_matches)
	cds_codons=cds_codons%>%GRangesList%>%unlist
	cds_codons$codon = names(cds_codons)%>%str_split('\\.')%>%map_chr(1)
	cds_codons
}

get_cov_rust_scores <- function(covgr,cds_codons,sampname){
	cdssampfpcov <- covgr%>%split(.,.$readlen)
	cdssampfpcov%>%map_df(.id='rl',function(cdsrlfpcov){
			rustcdsrlfpcov <- cdsrlfpcov
			rustcdsrlfpcov <- rustcdsrlfpcov>mean(rustcdsrlfpcov)
			nz_trs <- any(rustcdsrlfpcov)%>%names(.)[.]
			#
			cds_codons_nz = cds_codons%>%subset(seqnames%in%nz_trs)
			codtrs = seqnames(cds_codons_nz)%>%as.character
			cat('.')
			#calculate ro vals
			ro_cl = rustcdsrlfpcov[cds_codons_nz]%>%
				split(cds_codons_nz$codon)%>%
				lapply(as.matrix)%>%
				map(colMeans)
			#also get evals
			tr_rust_evals <- rustcdsrlfpcov%>%mean
			re_c <- tr_rust_evals[codtrs]%>%
				split(cds_codons_nz$codon)%>%
				map_dbl(mean)
			ro_cl <- ro_cl%>%map_df(.id='codon', enframe, 'position', 'ro_cl')
			re_c <- enframe(re_c, 'codon', 're_c')
			ro_cl%>%left_join(re_c, by='codon')
	})
}

get_sample_profs <- function(covgrs, cds_codons, n_wind_l_ext=45){
	stopifnot(!is.null(names(covgrs)))
	cdsfpcovlist <- lapply(covgrs,function(covgr) coverage(split(resize(covgr,1,'start'),covgr$readlen)))
	rust_roel <- mclapply(mc.cores=4,SIMPLIFY=F,cdsfpcovlist,cds_codons=cds_codons,get_cov_rust_scores)
	#https://www.nature.com/articles/ncomms12915#Sec10 see equation 3
	# of RUST paper
	frustprofilelist <- rust_roel%>%
	# frustprofilelist <- cn_norm%>%
		map_df(.id='set',.%>%map_df(.id='sample',.%>%bind_rows(.id='readlen')))
	#
	frustprofilelist %<>% mutate(position = position - 1 - (n_wind_l_ext))
		frustprofilelist %<>% filter(!codon %in% c('TAG','TAA','TGA'))
	frustprofilelist%<>%mutate(count = ro_cl/re_c)
	frustprofilelist$nreadlen <-frustprofilelist$readlen%>%as.numeric
	frustprofilelist$readlen%<>%str_replace('^(\\d)','rl\\1')
	frustprofilelist
}

#
get_kl_df <- function(frustprofilelist){
	frustprofilelist%>%
	filter(set=='all')%>%
	group_by(set,sample,nreadlen,position)%>%
	mutate(ro_cl = ro_cl/sum(ro_cl), re_c = re_c/sum(re_c))%>%
	summarise(KL=sum(ro_cl * log2(ro_cl/re_c)))
}

most_freq <- function(x) x%>% table%>%sort%>%names%>%as.numeric%>%tail(1)

select_offsets <- function(kl_df){
	if(method=='a_max'){
		kl_offsets <- kl_df%>%
			group_by(nreadlen,position)%>%
			summarise(sumKL = sum(KL))%>%
			filter(position< -6)%>%
			filter(position> -(nreadlen-6))%>%
			slice(which.max(sumKL))%>%
			mutate(p_offset = position+3)		
	}
	kl_offsets%>%select(nreadlen,position,p_offset)
}

plot_kl_dv <- function(kl_df, kl_offsets, selreadlens=NULL){
	stopifnot(c('position','nreadlen','KL','sample')%in%colnames(kl_df))
	stopifnot(c('nreadlen','sample','p_offset')%in%colnames(kl_offsets))
	#
	kl_offsets2plot <- kl_offsets%>%filter(nreadlen%in%selreadlens)
	#
	kl_df%>%
		filter(position< -3)%>%
		filter(position> -(nreadlen-6))%>%
		# separate(sample,c('fraction','genotype','rep'),remove=F)%>%
		filter(nreadlen%in%selreadlens)%>%
		{
			qplot(data=.,x=position,y=KL)+
			theme_bw()+
			facet_grid(nreadlen~sample)+
			scale_y_continuous('RUST KL divergence')+
			scale_x_continuous('5 read position relative to codon ')+
			geom_vline(data=kl_offsets2plot,aes(xintercept= p_offset-3),color=I('green'),linetype=2)+
			geom_vline(data=kl_offsets2plot,aes(xintercept= p_offset),color=I('blue'),linetype=2)+
			ggtitle("RUST KL divergence vs position")
		}
}

