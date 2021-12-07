
# /fast/AG_Ohler/dharnet/eif4f_pipeline/ext_data/gencode.v37.primary_assembly.annotation.gtf
system(str_interp("grep -e 'chr22[[:space:]]' ${gtf} > data/gcv37.anno.chr22.gtf"))


cov %>%
  subset(seqnames %in% chr22trs) %>%
  rtracklayer::export("test.bam")
cov %>%
  subset(seqnames %in% chr22trs) %>%
  rtracklayer::export("test.bam")
# covbak <- cov



mainchr22orfs <- anno$cdsgrl[anno$uORF %>%
  .[!.] %>%
  names()] %>%
  unlist() %>%
  subset(seqnames == "chr22")
uORFchr22orfs <- anno$cdsgrl[anno$uORF %>%
  .[.] %>%
  names()] %>%
  unlist() %>%
  subset(seqnames == "chr22")
chr22orfs <- c(mainchr22orfs, uORFchr22orfs) %>%
  names() %>%
  unique()
subanno <- subset_annotation(chr22orfs, anno)
codposdf <- get_codposdf(chr22orfs, anno, fafile)
usethis::use_data(codposdf)
codposdf %>% mutate(count = rbindom(rate))
chr22trs <- subanno$exonsgrl %>% names()

chr22_anno <- anno
usethis::use_data(chr22_anno)
rpfs <- get_readgr(testbam, chr22_anno)
usethis::use_data(rpfs)
usethis::use_data(offsets_df)

ms_df <- readr::read_tsv("/fast/AG_Ohler/dharnet/eif4f_pipeline/tables/sdr_ms_ribo_table.tsv")

ms_df <- ms_df %>% filter(gene_id %in% str_replace(gritpms$gene_id, "\\.\\d+", ""))

use_data(ms_df)
ah = AnnotationHub()
ah[['AH75191']]

data(rpfs)
trfstrands = strand(oldchr22$exonsgrl)%>%{.[as(rep(1,length(.)),'IntegerList')]}%>%unlist%>%as.character
names(trfstrands) = names(oldchr22$exonsgrl)
strand(rpfs) <- trfstrands[as.character(seqnames(rpfs))]
nrpfs <- rpfs%>%mapFromTranscripts(oldchr22$exonsgrl)
nrpfs <- nrpfs[width(nrpfs)==rpfs$readlen[nrpfs$xHits]]
nrpfs <- keepSeqlevels(nrpfs, seqlevels(chr22_anno$exonsgrl))
nrpfs <- unique(nrpfs)
nrpfs <- nrpfs%>%mapToTranscripts(chr22_anno$exonsgrl)
nrpfs <- nrpfs[width(nrpfs)==rpfs$readlen[nrpfs$xHits]]
tbam = 'inst/extdata/nchr22.bam'
nrpfs%>%rtracklayer::export(tbam)