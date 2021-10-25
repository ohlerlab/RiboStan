get_codposdf <- function(ribocovtrs, anno, fafile){
	#
	fafileob <- Rsamtools::FaFile(fafile)
	#
	cdsgrl <- anno%>%subset(type=='CDS')%>%
		split(.,.$transcript_id)
	#
	bestcdsseq = cdsgrl[ribocovtrs]%>%
		sort_grl_st%>%
		GenomicFeatures::extractTranscriptSeqs(x=fafileob,.)
	#
	codposdf = lapply(bestcdsseq,function(cdsseq){
		codonmat = codons(cdsseq)%>%
			{cbind(pos = .@ranges@start,as.data.frame(.))}%>%
			identity
	})
	codposdf%<>%bind_rows(.id='transcript_id')
	codposdf
}
get_psitecov <- function(sampled_cov_gr, offsets_df, anno){
	offsets <- sampled_cov_gr%>%
		# subset(readlen %in% offsets_df$readlen)%>%
		addphase(anno$cdsstarts)%>%
		{as.data.frame(mcols(.)[,c('readlen', 'phase')])}%>%
		left_join(
			offsets_df%>%select(readlen,phase,offset),
			by=c('readlen', 'phase')
		)%>%.$offset
	#
	psitecov <- sampled_cov_gr[!is.na(offsets)]%>%
		# subset(seqnames==zeropsitetr)
		resize(1, 'start')%>%
		shift(offsets%>%keep(Negate(is.na)))%>%
		coverage
	# psitecov%>%sum%>%min
	psitecov[anno$trspacecds[names(psitecov)]]
}
#
get_sitedf<-function(psitecov, anno, fafile,stop_codons=c('TAA','TAG','TGA')){
	#
	sitedf <- psitecov%>%
		lapply(.%>%matrix(nrow=3)%>%colSums)%>%
		stack%>%
		set_colnames(c('count','transcript_id'))
	#
	codposdf <- get_codposdf(names(psitecov), anno$anno, fafile)
	#
	stopifnot(
		all(unique(codposdf$transcript_id)==unique(sitedf$transcript_id))
	)
	#
	sitedf$p_codon<-codposdf$x
	sitedf$a_codon <- codposdf%>%split(.,.$transcript_id)%>%map(~lead(.$x))%>%
		unlist
	# psitecodons = codposdf%>%.[names(psitecov)]%>%bind_rows
	# #
	# sitedf<-sitedf%>%group_by(gene)%>%
	# 	mutate(phase=as_factor(((1:length(count))-1)%%3) )
	# sitedf$codon=NA
	# sitedf$codon[seq(1,nrow(sitedf),by=3)] = psitecodons$x
	# sitedf$codon[seq(2,nrow(sitedf),by=3)] = psitecodons$x
	# sitedf$codon[seq(3,nrow(sitedf),by=3)] = psitecodons$x
	# acod = psitecodons%>%
	# 	group_by(transcript_id)%>%mutate(acod=lead(x))%>%.$acod
	# sitedf$a_codon = NA
	# sitedf$a_codon[seq(1,nrow(sitedf),by=3)] = acod
	# sitedf$a_codon[seq(2,nrow(sitedf),by=3)] = acod
	# sitedf$a_codon[seq(3,nrow(sitedf),by=3)] = acod
	# sitedf <- sitedf%>%
		# filter(!codon %in% stop_codons,!a_codon %in% stop_codons)
	sitedf
}

get_codon_occs <- function(sampled_cov_gr, offsets_df,
		anno, n_genes=1000, method='linear'){
	#
	allpsitecov <- get_psitecov(sampled_cov_gr, offsets_df, anno)
	#
	#
	psitecov <- allpsitecov[allpsitecov%>%sum%>%`<`(10)]
	toptrs = psitecov%>%mean%>%sort%>%tail(n_genes)%>%names
	#
	psitecov <- psitecov[toptrs]
	#

	sitedf <- get_sitedf(psitecov, anno, fafile)
	library(MASS)
	#
	# message('faking data')
	# sitedf$count = 100
	# sitedf$count = sitedf$count * codonstrengths[sitedf$p_codon]
	# sitedf$count%<>%rpois(length(.),.)
	sitedf%<>%group_by(transcript_id)%>%mutate(trmean = mean(count))
	sitedf = sitedf%<>%filter(trmean!=0)
	
	sitedf

	#
	if(method =='linear'){
		nbfit = glm.nb(
			data=sitedf,
			formula= count ~ 0 + offset(log(trmean)) + p_codon + a_codon
		)
		#
		codon_occ_df <- broom::tidy(nbfit)%>%as.data.frame%>%
			mutate(codon=term%>%str_replace('[pa]_codon',''))%>%
			mutate(position=term%>%str_extract('^[pa]_codon'))%>%
			mutate(lower = estimate - 1.96*std.error)%>%
			mutate(upper = estimate + 1.96*std.error)
		codon_occ_df
		# codon_occ_df%>%filter(position=='p_codon')%>%
		# 	mutate(comp=codonstrengths[codon])%>%
		# 	{txtplot(.$comp,exp(.$estimate))}

	} else if(method=='RUST_glm'){
		sitedf%<>%mutate(rust=count > trmean,rtrmean = mean(rust))
		rustfit = glm(
			data=sitedf,
			formula= rust ~ 0 + offset(log(rtrmean)) + p_codon + a_codon,
			family='binomial'
		)
		rcodon_occ_df <- broom::tidy(rustfit)%>%as.data.frame%>%
			mutate(codon=term%>%str_replace('[pa]_codon',''))%>%
			mutate(position=term%>%str_extract('^[pa]_codon'))%>%
			mutate(lower = estimate - 1.96*std.error)%>%
			mutate(upper = estimate + 1.96*std.error)
		codon_occ_df = rcodon_occ_df
		# codon_occ_df%>%filter(position=='p_codon')%>%
		# 	mutate(comp=codonstrengths[codon])%>%
		# 	{txtplot::txtplot(.$comp,.$estimate)}
	} else if(method=='RUST'){
		sitedf%<>%mutate(rust=count > trmean,rtrmean = mean(rust))
			p_ests = sitedf%>%group_by(p_codon)%>%
			summarise(estimate = sum(rust)/sum(trmean))%>%
			setNames(c('codon','estimate'))%>%
			mutate(position = 'p_codon')
		a_ests = sitedf%>%group_by(p_codon)%>%
			summarise(estimate = sum(rust)/sum(trmean))%>%
			setNames(c('codon','estimate'))%>%
			mutate(position = 'a_codon')
		codon_occ_df = bind_rows(p_ests,a_ests)
		codon_occ_df%<>%mutate(upper=NA,lower=NA,p.value=NA)
		codon_occ_df
	}
	else{
		p_ests = sitedf%>%group_by(p_codon)%>%
			summarise(estimate = sum(count)/sum(trmean))%>%
			setNames(c('codon','estimate'))%>%
			mutate(position = 'p_codon')
		a_ests = sitedf%>%group_by(p_codon)%>%
			summarise(estimate = sum(count)/sum(trmean))%>%
			setNames(c('codon','estimate'))%>%
			mutate(position = 'a_codon')
		codon_occ_df = bind_rows(p_ests,a_ests)
		codon_occ_df%<>%mutate(upper=NA,lower=NA,p.value=NA)
		codon_occ_df
	}
	codon_occ_df%>%
		select(codon,position,estimate,upper,lower,p.value)
}



# #testing various footprint methods - they all line up pretty well
# cor.test(codon_occ_df$estimate,rcodon_occ_df$estimate)

# linear_codon_occ_df <- get_codon_occs(sampled_cov_gr, offsets_df, anno,
# 	n_genes=1000, method='linear')
# rust_codon_occ_df <- get_codon_occs(sampled_cov_gr, offsets_df, anno,
# 	n_genes=1000, method='RUST_glm')
# std_codon_occ_df <- get_codon_occs(sampled_cov_gr, offsets_df, anno,
# 	n_genes=1000, method='standard')
# srust_codon_occ_df <- get_codon_occs(sampled_cov_gr, offsets_df, anno,
# 	n_genes=1000, method='RUST')
# # std_codon_occ_df$position%<>%paste0(.,'_codon')
# rust_codon_occ_df%>%filter(!codon%in%c('TAG','TAA','TGA'))%>%filter(position=='a_codon')%>%.$estimate%>%txtdensity
# rust_codon_occ_df%>%filter(!codon%in%c('TAG','TAA','TGA'))%>%filter(position=='p_codon')%>%.$estimate%>%txtdensity

# cor.test(linear_codon_occ_df$estimate,rust_codon_occ_df$estimate)
# txtplot(linear_codon_occ_df$estimate,rust_codon_occ_df$estimate)

# linear_codon_occ_df%>%
# 	inner_join(rust_codon_occ_df,by=c('codon','position'))%>%
# 	filter(position=='p_codon')%>%
# 	filter(!codon%in%c('TAG','TAA','TGA'))%>%
# 	{txtplot(.$estimate.x,.$estimate.y)}


# linear_codon_occ_df%>%
# 	inner_join(std_codon_occ_df,by=c('codon','position'))%>%
# 	filter(position=='p_codon')%>%
# 	filter(!codon%in%c('TAG','TAA','TGA'))%>%
# 	{txtplot(.$estimate.x,.$estimate.y)}

# linear_codon_occ_df%>%
# 	inner_join(srust_codon_occ_df,by=c('codon','position'))%>%
# 	filter(position=='p_codon')%>%
# 	filter(!codon%in%c('TAG','TAA','TGA'))%>%
# 	{txtplot(.$estimate.x,.$estimate.y)}

# rust_codon_occ_df%>%
# 	inner_join(srust_codon_occ_df,by=c('codon','position'))%>%
# 	filter(position=='p_codon')%>%
# 	filter(!codon%in%c('TAG','TAA','TGA'))%>%
# 	{txtplot(.$estimate.x,.$estimate.y)}
#this probably needs to take an actual model object.
get_predicted_codon_occs <- function(anno, fafile, codon_model, ribocovtrs=NULL){
	if(is.null(ribocovtrs))ribocovtrs <- anno$transcript_id%>%unique
	#
	cdsgrl <- anno$anno%>%subset(type=='CDS')%>%
		split(.,.$transcript_id)
	#
	fafileob <- Rsamtools::FaFile(fafile)
	#
	elongtrs <- unique(anno$anno$transcript_id)
	trsplit = elongtrs%>%split(floor(seq_along(.)/1000))
	tr_elong = map_df(trsplit,function(elongtrs){
		#
		bestcdsseq = cdsgrl[elongtrs]%>%
			sort_grl_st%>%
			GenomicFeatures::extractTranscriptSeqs(x=fafileob,.)

		codonfreqdf = bestcdsseq%>%oligonucleotideFrequency(3,step = 3)%>%
			set_rownames(names(bestcdsseq))%>%
			{.=./rowSums(.);.}%>%
			as.data.frame%>%rownames_to_column('transcript_id')%>%
			pivot_longer(-transcript_id,names_to='codon',values_to='freq')

		ipos = 'a_codon'
		positions = unique(codon_model$position)
		poselongs <- lapply(positions,function(ipos){
			poselong = codonfreqdf%>%left_join(
				codon_model%>%filter(position==ipos)%>%
				select(codon,estimate),
				by='codon'
			)%>%
			group_by(transcript_id)%>%
			summarise(estimate=weighted.mean(estimate,freq,na.rm=T))
			poselong
		})
		tr_elong = poselongs%>%purrr::reduce(left_join,by='transcript_id')%>%
			{tibble(transcript_id=.[[1]],sum_occ = rowSums(.[,-1]))}
		tr_elong
	})
	tr_elong
}


gene_level_elong<-function(tr_elong, tpms, anno){
	trgiddf = anno$anno%>%mcols%>%
		.[,c('transcript_id','gene_id')]%>%
		as.data.frame%>%
		distinct
		#
	gn_elong <- tr_elong%>%
			left_join(trgiddf,by='transcript_id')%>%
			group_by(gene_id)%>%
			left_join(enframe(tpms, 'transcript_id', 'tpm'),by='transcript_id')%>%
			ungroup%>%
			mutate(tpm = replace_na(tpm,min(tpm,na.rm=T)*0.01))%>%
			group_by(gene_id)%>%
			summarise(sum_occ = weighted.mean(sum_occ, tpm, na.rm=T))
	gn_elong%>%select(gene_id,sum_occ)
}

#
# probgids = gn_elong%>%filter(is.na(sum_occ))%>%pluck('gene_id')
# okay so some of them have 0 expression gene wise
# gtpms%>%filter(gene_id%in%probgids)
# gns_w_expr = gtpms%>%filter(expr!=0)%>%pluck('gene_id')
#
# probgids_w_expr = gn_elong%>%filter(is.na(sum_occ),gene_id%in%gns_w_expr)%>%pluck('gene_id')

