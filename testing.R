
{
knitr::opts_chunk$set(root.dir = here::here(),eval=TRUE,cache=FALSE,echo=FALSE,warning = FALSE,message = FALSE,include=TRUE,
                      fig.width =7,fig.height=7,out.width=700,out.height=700,dev='png')
library(dplyr)
library(stringr)
library(ggplot2)
library(GenomicRanges)
library(ORFik)
shift<-GenomicRanges::shift
#' @importFrom tidyr replace_na
replace_na <- tidyr::replace_na

dpss <- multitaper::dpss
spec.mtm <- multitaper::spec.mtm
dropFreqs <- multitaper::dropFreqs
filter <- dplyr::filter
mutate <- dplyr::mutate
select <- dplyr::select
pivot_longer <- tidyr::pivot_longer
rownames_to_column <- tibble::rownames_to_column
}
  #
gtf <- '../cortexomics/ext_data/gencode.vM12.annotation.gtf'
ribobam <- '../cortexomics/pipeline/star_transcript/data/E13_ribo_1/E13_ribo_1.bam'

gtf <- '../eif4f_pipeline/ext_data/gencode.v37.primary_assembly.annotation.gtf'
ribobam <- '../eif4f_pipeline/pipeline/merged_pc_wt.sort.bam'
source('R/anno.R')
source('R/get_tpms.R')
source('R/get_offsets.R')
source('R/find_ORFs.R')
#
#
ribobams <- Sys.glob('/fast/AG_Ohler/dharnet/eif4f_pipeline/pipeline/star/pc/data/4E_-IAA_rep*/*.bam')
ribobam <- ribobams[1]
fafile = '/fast/AG_Ohler/dharnet/eif4f_pipeline/GRCh38.primary_assembly.genome.fa'
# make_ext_fasta(gtf=gtf,fa=fafile,out="foo.fa")
#
anno <- rtracklayer::import(gtf)
anno <- filter_anno(anno, fafile)

txdb <- GenomicFeatures::makeTxDbFromGFF(gtf)

library(GenomicFeatures)

fiveutrs = fiveUTRsByTranscript(txdb, use.names=T)
alluORFs <- findUORFs(fiveutrs, fafile)


{
#
source('R/anno.R')
fafile = '/fast/AG_Ohler/dharnet/eif4f_pipeline/GRCh38.primary_assembly.genome.fa'
anno <- load_annotation(gtf, fafile)
}

source('R/get_offsets.R')
source('R/get_ritpms.R')
#
cov <- get_readgr(ribobam, anno)
offsets_df <- get_offsets(cov, anno)
psites <- get_psite_gr(get_psite_gr, offsets_df, anno)
#
{
source('R/get_ritpms.R')
#
ritpms <- get_ritpms(ritpm_opt)
gritpms = gene_level_expr(ritpms, anno)
}

{
source(('R/model_TE.R'))
source('R/get_offsets.R')
source('R/get_ritpms.R')
ribocovtrs <- ritpms%>%Filter(f=function(x)x>10)%>%names
mainORFs <- anno$trgiddf%>%filter(!uORF)%>%.$transcript_id
ribocovtrs <- intersect(ribocovtrs, mainORFs)
#save.image('dev.Rdata')
#
rust_codon_occ_df <- get_codon_occs(psites, offsets_df, anno,
  n_genes=1000, method='RUST_glm')
tr_elong = get_predicted_codon_occs(anno, fafile, rust_codon_occ_df)
gn_elong = gene_level_elong(tr_elong, ritpms, anno)
}

{
selreadlens=27:33
source('R/kl_div.R')
cds_codons <- get_cds_codons(anno)
covgrs = list(sample1=sampled_cov_gr)
fprustprofilelist <- get_sample_profs(covgrs, cds_codons, n_wind_l_ext=45, )
kl_df<-get_kl_df(fprustprofilelist)
kl_offsets <- select_offsets(kl_df)
kl_offsets%>%select(p_offset,length=nreadlen)%>%readr::write_tsv('offsets_rustvar.tsv')
#
kl_div_plot <- plot_kl_dv(kl_df, kl_offsets, selreadlens)
#
# pdf<-grDevices::pdf
# dir.create('plots')
# plotfile='plots/rust_fppos_vs_codon_variance.pdf'
# pdf(plotfile,w=12,h=1+3*length(selreadlens))
kl_div_plot%>%print
# dev.off()
# message(normalizePath(plotfile))
#
allcodondt <- export_codon_dts(fprustprofilelist, kl_offsets)
allcodondt %>% readr::write_tsv('allcodondt.tsv')
#
}

{

}

#annotation
#data loading
#offsets
#Metaplots
#Kl plot
#quantification
#Detect Periodicity
#Gene level Elong
#Comparison to MS data

#Now compare calculated elongation rates to the proteomic derived values
sdr_df <- readr::read_tsv('/fast/AG_Ohler/dharnet/eif4f_pipeline/tables/sdr_ms_ribo_table.tsv')
#
library(txtplot)
sdr_df%>%left_join(gritpms%>%mutate_at('gene_id',str_replace,'\\.\\d+$',''))%>%
  filter(is.finite(log2(expr)))%>%
  filter(ribo>0)%>%
  # {txtplot(.$ribo,log2(.$expr))}%>%
  {cor.test(.$ribo,log2(.$expr))}

sdr_df%>%left_join(gritpms%>%mutate_at('gene_id',str_replace,'\\.\\d+$',''))%>%
  filter(is.finite(log2(expr)))%>%
  # filter(expr>0)%>%
  # {txtplot(.$ribo,log2(.$expr))}%>%
  {cor.test(.$ms,log2(.$expr))}

sdr_df%>%left_join(gritpms%>%mutate_at('gene_id',str_replace,'\\.\\d+$',''))%>%
  filter(is.finite(log2(expr)))%>%
  # filter(expr>0)%>%
  # {txtplot(.$ribo,log2(.$expr))}%>%
  {cor.test(.$ms,log2(.$ribo))}


sdr_df%>%left_join(gn_elong%>%mutate_at('gene_id',str_replace,'\\.\\d+$',''))%>%
  filter(is.finite(log2(sum_occ)))%T>%
  {txtplot(.$sdr,log2(.$sum_occ))}%>%
  {cor.test(.$sdr,log2(.$sum_occ))}

#does the ribo correlation with my 'ritpms' here?
sdr_df%>%left_join(enframe(ritpms,'gene_id','ritpm'))%>%{cor.test(.$ribo,.$ritpm)}
sdr_df%>%left_join(enframe(ritpms,'gene_id','ritpm'))%>%{cor.test(.$ribo,.$ritpm)}


rust_profiles <- sampled_cov_gr%>%create_rust_profiles
#
elong_rates <- rust_profiles%>%get_elong_rates()