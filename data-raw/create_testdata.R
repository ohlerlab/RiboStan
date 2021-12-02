
/fast/AG_Ohler/dharnet/eif4f_pipeline/ext_data/gencode.v37.primary_assembly.annotation.gtf
system(str_interp("grep -e 'chr22[[:space:]]' ${gtf} > data/gcv37.anno.chr22.gtf"))


cov%>%subset(seqnames%in%chr22trs)%>%rtracklayer::export('test.bam')
cov%>%subset(seqnames%in%chr22trs)%>%rtracklayer::export('test.bam')
# covbak <- cov



mainchr22orfs <- anno$cdsgrl[anno$uORF%>%.[!.]%>%names]%>%unlist%>%subset(seqnames=='chr22')
uORFchr22orfs <- anno$cdsgrl[anno$uORF%>%.[.]%>%names]%>%unlist%>%subset(seqnames=='chr22')
chr22orfs<-c(mainchr22orfs,uORFchr22orfs)%>%names%>%unique
subanno <- subset_annotation(chr22orfs,anno)
codposdf <- get_codposdf(chr22orfs, anno, fafile)
usethis::use_data(codposdf)
codposdf %>% mutate(count = rbindom(rate))
chr22trs <- subanno$exonsgrl%>%names

chr22_anno <- anno
usethis::use_data(chr22_anno)
rpfs <- get_readgr(testbam, chr22_anno)
usethis::use_data(rpfs)
usethis::use_data(offsets_df)

ms_df <- readr::read_tsv('/fast/AG_Ohler/dharnet/eif4f_pipeline/tables/sdr_ms_ribo_table.tsv')

ms_df= ms_df%>%filter(gene_id %in% str_replace(gritpms$gene_id,'\\.\\d+',''))

use_data(ms_df)