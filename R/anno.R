#now only those which have M at the start and '*' at the end
hasMstart <- function(cdsgrl, fafileob){
	cdsseqstarts = cdsgrl%>%
		sort_grl_st%>%
		resize_grl(3,'start')%>% 
		GenomicFeatures::extractTranscriptSeqs(x=fafileob,.)%>%
		Biostrings::translate(.)
	cdsseqstarts[ribocovtrs]=='M'
}

get_trspace_cds <-function(cdsgrl,filt_anno){
	exonsgrl <- filt_anno%>%
		subset(type=='exon')%>%
		split(.,.$transcript_id)%>%
		.[names(cdsgrl)]
	#now lift cds to exons space
	trspacecds <- GenomicFeatures::pmapToTranscripts(
		cdsgrl,
		exonsgrl[names(cdsgrl)])
	trspacecds <- trspacecds%>%unlist
	names(trspacecds)==names(cdsgrl)
	trspacecds
}

#get the grl 
get_cdsgrl <- function(filt_anno, ignore_orf_validity){
	#find which cds are multiples of 3bp
	cdsgrl <- filt_anno%>%subset(type=='CDS')%>%split(.,.$transcript_id)
	is3bp = cdsgrl%>%width%>%sum%>%`%%`(3)%>%`==`(0)	
	cdsgrl <- cdsgrl[is3bp]
	message(str_interp('filtered out ${sum(!is3bp)} ORFs for not being multiples of 3bp long'))
	#chceck if the cds includes the stop codon
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
	#now determine if the annotations 'cds' include stop codons
	#if they do, fix that.
	end_stop = table(subseq(cdsseqends,1,1))%>%sort%>%
	{./sum(.)}%>%.['*']%>%`>`(0.5)
	if(is.na(end_stop)) end_stop = FALSE
	end_plusone_stop = table(subseq(cdsseqends,2,2))%>%sort%>%
		{./sum(.)}%>%.['*']%>%`>`(0.5)
	stopifnot(end_stop|end_plusone_stop)
	if(end_stop) cdsgrl <- cdsgrl%>%resize_grl(sum(width(.))-3,'start')
	#
	endseq = if(end_plusone_stop){ subseq(cdsseqends,2,2) }else{
	 subseq(cdsseqends,1,1)}
	hasstop = endseq[ribocovtrs]=='*'
	if(!ignore_orf_validity){
		cdsgrl <- cdsgrl[hasstop]
		message(str_interp('filtered out ${sum(!hasstop)} ORFs for not ending with *'))
	}
	hasM <- hasMstart(cdsgrl)
	if(!ignore_orf_validity){
		cdsgrl <- cdsgrl[hasM]
		message(str_interp('filtered out ${sum(!hasstop)} ORFs for not starting with M'))
	}
	message(str_interp('${length(cdsgrl)} ORFs left'))
	cdsgrl

}
################################################################################
########## 
################################################################################
filter_anno <- function(anno, fafile, ignore_orf_validity=F){
	filt_anno <- anno
	fafileob = Rsamtools::FaFile(fafile)
	require(Biostrings)

	#get the cds not including stop codons, possibly filtering for valid orfs
	cdsgrl<-get_cdsgrl(filt_anno, fafileob, ignore_orf_validity)
	#
	#subset cds and anno with these
	filt_anno = filt_anno%>%
		subset(type!='CDS')%>%
		subset(transcript_id%in%names(cdsgrl))
	filt_anno = c(filt_anno,unlist(cdsgrl))
	#
	ribocovtrs = names(cdsgrl)
	#
	trspacecds <- get_trspace_cds(cdsgrl, exonsgrl)
	#
	cdsstarts = trspacecds%>%start%>%setNames(names(trspacecds))
	#
	exonsgrl <- filt_anno%>%
		subset(type=='exon')%>%
		split(.,.$transcript_id)%>%
		.[names(cdsgrl)]
	outanno = list(
		ribocovtrs=ribocovtrs,
		trspacecds=trspacecds,
		cdsgrl=cdsgrl,
		exonsgrl= exonsgrl)
	outanno = c(outanno,
		list(cdsstarts = outanno$trspacecds%>%start%>%
			setNames(names(outanno$trspacecds)),
		cds_prestop_st = outanno$trspacecds%>%end%>%`-`(2)%>%
			setNames(names(outanno$trspacecds)),
		anno = filt_anno%>%subset(transcript_id%in%ribocovtrs)		
	))
	return(outanno)
}

make_ext_fasta <- function(gtf, fasta, outfasta){
	stopifnot({cat('testing',file=outfasta);file.remove(outfasta)})

	stopifnot(gtf%>%str_detect('\\.(gtf)$'))
	stopifnot(gtf%>%file.exists)

	stopifnot(outfasta%>%str_detect('\\.(fasta|fa)$'))
	outprefix = outfasta%>%str_replace('\\.(fasta|fa)$','')

	#get our filtered annotation
	anno <- filter_anno(rtracklayer::import(gtf), fasta)
	exonsgrl = anno$exonsgrl
	cdsexonsgrl = anno$cdsgrl

	#get exons for our cds
	cdsexonsgrl%<>%sort_grl_st
	#get an object representing the CDS In transript space
	cdstrspace = pmapToTranscripts(filtcdsgrl,cdsexonsgrl[fmcols(filtcdsgrl,transcript_id)])
	#ensure all cds map cleanly to the exons
	stopifnot(cdstrspace%>%elementNROWS%>%`==`(1))
	cdsstartpos = setNames(start(cdstrspace@unlistData),names(cdstrspace))
	endpos = sum(width(cdsexonsgrl))-end(cdstrspace@unlistData)
	#expand our first exon when needed
	startposexpansion = pmax(0,fpext - cdsstartpos + 1)
	#expand/trim the 5' end of the exons
	cdsexonsgrl@unlistData[start(cdsexonsgrl@partitioning)]%<>%resize(width(.)+startposexpansion,'end')
	#expand or trim the last exon when needed
	endposexpansion = pmax(0,tpext - endpos)
	cdsexonsgrl@unlistData[cdsexonsgrl@partitioning@end]%<>%resize(width(.)+endposexpansion,'start')

	{
	#now map our cds to that
	cds_exptrspc = pmapToTranscripts(filtcdsgrl,cdsexonsgrl)
	stopifnot(cds_exptrspc%>%elementNROWS%>%`==`(1))

	expcds_exptrspc=cds_exptrspc
	stopifnot(!any(expcds_exptrspc%>%elementNROWS%>%`>`(1)))
	expcds_exptrspc%<>%unlist
	#expand our cds exons
	expcds_exptrspc%<>%resize(width(.)+fpext,'end',ignore.strand=TRUE)
	#and expand the 3' ends
	expcds_exptrspc%<>%resize(width(.)+tpext,'start',ignore.strand=TRUE)
	#now back to genome space
	expcdsgenspace = spl_mapFromTranscripts(expcds_exptrspc,cdsexonsgrl)
	expcdsgenspace = split(expcdsgenspace,names(expcdsgenspace))
	#get the sequences
	isoutofbds = any(is_out_of_bounds(expcdsgenspace,seqinfo(fafile)))
	message(str_interp('Excluded ${sum(isoutofbds)} genes because they extended beyond chromosomal boundaries'))
	expcdsgenspace <- expcdsgenspace[!isoutofbds]
	cdsexonsgrl <- cdsexonsgrl[names(expcdsgenspace)]
	cds_exptrspc <- cds_exptrspc[names(expcdsgenspace)]
	expcdsgenspaceseq <- 
		expcdsgenspace%>%
		sort_grl_st%>%
		extractTranscriptSeqs(x=fafile)

	}
	cdslens = sum(width(cdsgrl))[names(expcdsgenspace)]
	{
	fastanames <- paste(sep='|',
		fmcols(cdsexonsgrl,transcript_id),
		fmcols(cdsexonsgrl,gene_id),
		fmcols(cdsexonsgrl,havana_gene),
		fmcols(cdsexonsgrl,havana_transcript),
		fmcols(cdsexonsgrl,transcript_name),
		fmcols(cdsexonsgrl,gene_name),
		sum(width(cdsexonsgrl)),
		paste0('UTR5:1-',fpext),
		paste0('CDS:',fpext+1,'-',fpext+cdslens),
		paste0('UTR3:',1+fpext+cdslens,'-',sum(width(expcdsgenspace))),
		'|')

	names(expcdsgenspaceseq) <- fastanames

	trcdscoordsfile <- paste0(outprefix,'_trcds.tsv')
	as.data.frame(unlist(cds_exptrspc))%>%select(seqnames,start,end)%>%write_tsv(trcdscoordsfile)
	message(normalizePath(trcdscoordsfile,mustWork=TRUE))
	}

	{
		#also write our cds coordinates to disk in the new trspace
		new_trspc_anno <- c(
			GRanges(as.character(seqnames(cds_exptrspc)),IRanges(1,sum(width(cdsexonsgrl))))%>%{.$type='exon';.},
			cds_exptrspc%>%unlist%>%{.$type='CDS';.}
		)
		new_trspc_anno%>%export(paste0(outprefix,'_trspaceanno.gtf'))

		#write the expanded cds exon sequences to disk
		writeXStringSet(expcdsgenspaceseq,outfasta)
		message(normalizePath(outfasta),mustWork=TRUE)

		#now make fasta file with shorter transcript names
		shortnamefile = paste0(outprefix,'.shortheader.fa')
		system(str_interp("sed -e 's/|.*$//' ${outfasta} > ${shortnamefile}"))
		message(normalizePath(shortheaderfasta,mustWork=TRUE))

	}
	return(outfasta)
}