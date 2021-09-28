


if(!file.exists('data/txcodposdf.rds')) {
	fafileob <- Rsamtools::FaFile(fafile)
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
	codposdf%>%saveRDS(here('data/txcodposdf.rds'))
}else{
	codposdf<-readRDS(here('data/txcodposdf.rds'))
}

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
		addphase(cdsstarts)%>%
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

psitecov <- get_psitecov(sampled_cov_gr, offsets_df, anno)
psitecov <- psitecov[psitecov%>%sum%>%`<`(10)]

#debug
zeropsitetrs <- psitecov%>%sum%>%keep(equals,0)%>%names
zeropsitetr=zeropsitetrs[1]
tpms[zeropsitetr]
cov%>%subset(seqnames==zeropsitetr)
trspacecds[zeropsitetr]
sampled_cov_gr%>%subset(seqnames==zeropsitetr)
tpms[zeropsitetr]


################################################################################
########## 
################################################################################

testr='ENST00000598261.2'



psitecov <- psitecov%>%sample(100)

sitedf <- get_sitedf(psitecov%>%sample(100), anno, fafile)
library(MASS)

codglmfit = glm.nb(
	data=sitedf,
	formula= count ~ 0 + transcript_id + p_codon+a_codon)

glmfits = lapply(sections,function(section){
	# 
	tnt_entrop_cds <- sum(psitecov)[is3nt]%>%order(decreasing=TRUE)%>%.[section]
	#
	testpsitecovs = c(psitecov[is3nt][tnt_entrop_cds])
	#
	sitedf = get_sitedf(testpsitecovs,codposdf%>%split(.,.$transcript_id))
	#
	library(MASS)
	message('fit')
	codglmfit = glm.nb(
		data=sitedf,
		formula= count ~ 0 + gene + codon+a_codon+phase)
	message('done')
	codglmfit
})

