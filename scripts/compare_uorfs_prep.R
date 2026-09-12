#!/usr/bin/env Rscript
# Prepare Ribostan translation-supported uORFs + genomic BAM for RiboTaper comparison.
#
# Compares periodicity-filtered Ribostan uORFs to RiboTaper periodic ORFs only.
# Raw sequence-scan candidates are exported separately for diagnostics.
args <- commandArgs(trailingOnly = TRUE)
workdir <- if (length(args) >= 1) args[[1]] else "."
setwd(workdir)

suppressPackageStartupMessages({
  library(Ribostan)
  library(AnnotationHub)
  library(Biostrings)
  library(rtracklayer)
  library(GenomicRanges)
  library(GenomicFeatures)
  library(Rsamtools)
  library(magrittr)
  library(dplyr)
})

anno_file <- file.path(workdir, "test.gc32.gtf")
fafile    <- file.path(workdir, "chr22.fa")
if (!file.exists(anno_file)) {
  ah <- AnnotationHub()
  gencode32 <- ah[["AH75191"]]
  seqlevels(gencode32) <- "chr22"
  export(gencode32, anno_file)
}
if (!file.exists(fafile)) {
  fa_gz <- file.path(workdir, "chr22.fa.gz")
  ucsc_url <- "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr22.fa.gz"
  if (!file.exists(fa_gz)) {
    message("Downloading chr22.fa.gz from UCSC...")
    utils::download.file(ucsc_url, fa_gz, mode = "wb", quiet = TRUE)
  }
  seq <- Biostrings::readDNAStringSet(fa_gz)
  writeXStringSet(seq, fafile)
}

message("load_annotation(add_uorfs=TRUE)...")
chr22_anno <- load_annotation(anno_file, fafile, add_uorfs = TRUE)

uorf_to_df <- function(anno, label = "all") {
  uorf_ids <- names(anno$uORF)[anno$uORF]
  cds_grl  <- anno$cdsgrl[uorf_ids]
  tr_orfs  <- anno$trspacecds[uorf_ids]
  qvals <- mcols(tr_orfs)$q.value
  pvals <- mcols(tr_orfs)$p.value
  data.frame(
    set = label,
    orf_id = names(cds_grl),
    transcript_id = sub("_[0-9]+$", "", names(cds_grl)),
    seqnames = vapply(cds_grl, function(g) as.character(GenomicRanges::seqnames(g)[1]), character(1)),
    start = vapply(cds_grl, function(g) min(GenomicRanges::start(g)), integer(1)),
    end = vapply(cds_grl, function(g) max(GenomicRanges::end(g)), integer(1)),
    strand = vapply(cds_grl, function(g) as.character(GenomicRanges::strand(g)[1]), character(1)),
    width_nt = as.integer(GenomicRanges::width(tr_orfs)),
    width_codons = (as.integer(GenomicRanges::width(tr_orfs)) - 1L) / 3L,
    start_codon = as.character(unlist(S4Vectors::mcols(tr_orfs)$start_codon)),
    overlaps_cds = as.logical(unlist(S4Vectors::mcols(tr_orfs)$overlaps_cds)),
    p_value = as.numeric(pvals),
    q_value = as.numeric(qvals),
    n_psites = as.integer(unlist(S4Vectors::mcols(tr_orfs)$n_psites)),
    stringsAsFactors = FALSE
  )
}

candidate_df <- uorf_to_df(chr22_anno, "sequence_candidate")
write.table(candidate_df, file.path(workdir, "ribostan_uorfs_candidates.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
message("Ribostan sequence candidates: ", nrow(candidate_df))

testbam <- system.file("extdata", "nchr22.bam", package = "Ribostan", mustWork = TRUE)
rpfs <- get_readgr(testbam, chr22_anno)
message("RPFs: ", length(rpfs))

readlens <- sort(unique(as.integer(rpfs$readlen)))
offsets_df <- data.frame(
  readlen = readlens,
  p_offset = rep(12L, length(readlens)),
  phase = rep(0L, length(readlens)),
  stringsAsFactors = FALSE
)
message("Using fixed p_offset=12 for read lengths: ", paste(readlens, collapse = ", "))
writeLines(as.character(offsets_df$readlen), file.path(workdir, "ribo_readlengths.txt"))
writeLines(as.character(offsets_df$p_offset), file.path(workdir, "ribo_offsets.txt"))

psites <- get_psite_gr(rpfs, offsets_df, chr22_anno)
message("P-sites: ", length(psites))

filt_anno <- periodicity_filter_uORFs(psites, chr22_anno, remove = TRUE, n_cores = 1)
periodic_df <- uorf_to_df(filt_anno, "periodic")
periodic_df$set <- "periodic"

saveRDS(periodic_df, file.path(workdir, "ribostan_uorfs.rds"))
write.table(periodic_df, file.path(workdir, "ribostan_uorfs.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
message("Ribostan translation-supported uORFs (q < 0.05, n_psites >= 6): ", nrow(periodic_df))

tx_ids <- as.character(seqnames(rpfs))
exons  <- chr22_anno$exonsgrl[tx_ids]
mapped <- GenomicFeatures::pmapFromTranscripts(rpfs, exons)
if ("hit" %in% colnames(S4Vectors::mcols(mapped@unlistData))) {
  hit <- mapped@unlistData$hit
  mapped <- mapped[IRanges::relist(hit, mapped)]
  mapped@unlistData$hit <- NULL
}
gen_rpfs <- unlist(mapped, use.names = FALSE)
S4Vectors::mcols(gen_rpfs) <- S4Vectors::mcols(rpfs)
names(gen_rpfs) <- names(rpfs)
seqlevels(gen_rpfs) <- seqlevels(chr22_anno$cdsgrl)

gen_bam <- file.path(workdir, "nchr22_genomic.bam")
export(gen_rpfs, gen_bam)
indexFa(fafile)
message("prep done")
