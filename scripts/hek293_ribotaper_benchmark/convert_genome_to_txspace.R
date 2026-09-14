#!/usr/bin/env Rscript
# Project a genome-aligned, primary-only-filtered BAM (RiboTaper-style
# alignment; see align_genome_ribotaper_style.sh) into the transcript-space
# read GRanges Ribostan's periodicity pipeline expects, using the CCDS/APPRIS
# annotation's exon structure (anno$exonsgrl) as the projection target.
#
# A read's single genomic primary locus can legitimately land inside more
# than one transcript's exons when isoforms share the same exon (e.g. a
# shared 5' UTR) - mapToTranscripts() emits one output row per compatible
# transcript in that case. That is intentional here: it mirrors Ribostan's
# own multimap="all" transcript-space behaviour, so the only thing this
# genome-first route changes relative to the transcriptome-BAM route is
# removing *true* genomic multi-mappers (repeats/paralogs) before that
# broadcast, via STAR's own alignment + primary-selection logic, rather than
# via Ribostan's NH-tag-based multimap filtering.
args <- commandArgs(trailingOnly = TRUE)
workdir <- if (length(args) >= 1) args[[1]] else "."
setwd(workdir)

args_all <- commandArgs(trailingOnly = FALSE)
script_path <- sub("--file=", "", args_all[grep("--file=", args_all)])
bench_dir <- if (length(script_path)) dirname(script_path) else "."
source(file.path(bench_dir, "hek293_gtf_config.R"), local = TRUE)

suppressPackageStartupMessages({
  library(Ribostan)
  library(GenomicAlignments)
  library(GenomicFeatures)
  library(GenomicRanges)
  library(Rsamtools)
})

gtf_cfg <- resolve_hek293_gtf(workdir)
gtf <- gtf_cfg$path
fafile <- file.path(workdir, "GRCh37.p13.genome.fa")
bam <- file.path(workdir, Sys.getenv("HEK293_GENOME_BAM",
                                     "hek293_genome_rrnafilt_primary.bam"))
out_rds <- file.path(workdir, Sys.getenv("HEK293_PROJECTED_RDS",
                                         "rpfs_genome_projected_rrnafilt_matched.rds"))

stopifnot(file.exists(gtf), file.exists(fafile), file.exists(bam))

message("Loading annotation (", gtf_cfg$name, ") ...")
# Projection only needs transcript exon structures. Discovering every uORF
# candidate here is expensive and cannot change the projected read coordinates.
anno <- load_annotation(gtf, fafile, add_uorfs = FALSE)

message("Reading primary-only genome BAM: ", bam)
flags <- Rsamtools::scanBamFlag(isUnmappedQuery = FALSE)
# Keep spliced (N-in-CIGAR) primary alignments. RiboTaper's genome pipeline
# does not drop these; simpleCigar=TRUE (as in read_ribobam) removes ~14% of
# genome primary reads and was a remaining data mismatch vs RiboTaper.
bparam <- Rsamtools::ScanBamParam(flag = flags, tag = "NH")
gal <- GenomicAlignments::readGAlignments(bam, use.names = TRUE, param = bparam)
message("Genome-aligned primary reads: ", length(gal))

get_readlens <- function(readlens) {
  # Mirrors Ribostan's own (unexported) get_readlens(): keep the middle 95%
  # of the length distribution by cumulative fraction, trimming both tails -
  # not just the bottom tail - so this matches the same read-length window
  # the transcriptome-BAM comparison run used.
  tab <- table(as.numeric(readlens))
  frac <- cumsum(tab) / sum(tab)
  keep <- frac[frac > 0.025]
  keep <- keep[(1 - keep) > 0.025]
  as.integer(names(keep))
}
readlens <- get_readlens(GenomicAlignments::qwidth(gal))
wfilt <- GenomicAlignments::qwidth(gal) <= max(readlens) &
  GenomicAlignments::qwidth(gal) >= min(readlens)
gal <- gal[wfilt]
message("After read-length filter [", min(readlens), ",", max(readlens),
        "]: ", length(gal), " reads")

readlen <- GenomicAlignments::qwidth(gal)
read_names <- names(gal)
gr <- as(gal, "GRanges")
seqlevels(gr) <- sub("\\|.*", "", seqlevels(gr))

# Keep projection targets transcript-unique defensively.
exonsgrl_uniq <- anno$exonsgrl[!duplicated(names(anno$exonsgrl))]
message("Projecting into transcript space against ", length(exonsgrl_uniq),
        " annotated transcripts ...")
mapped <- GenomicFeatures::mapToTranscripts(gr, exonsgrl_uniq, ignore.strand = FALSE)
message("Projected alignments (read x compatible-transcript rows): ", length(mapped))
message("Unique reads represented: ", length(unique(mapped$xHits)), " / ",
        length(gr), " input primary reads (",
        sprintf("%.1f", 100 * length(unique(mapped$xHits)) / length(gr)),
        "% land inside an annotated exon of a kept transcript)")

nh <- S4Vectors::mcols(gal)$NH
nh <- if (is.null(nh)) rep(1L, length(gal)) else as.integer(ifelse(is.na(nh), 1L, nh))
mcols(mapped)$readlen <- readlen[mapped$xHits]
# Keep the *genomic* NH (how many genomic loci STAR found for this read
# before -F 0x100 kept only the primary one) rather than forcing NH=1.
# This is the number RiboTaper's own multi-mapping/coverage QC
# (pct_covered_onlymulti/pct_covered_all) is calibrated against - true
# repeat/paralog ambiguity - as opposed to the isoform-sharing "multimapping"
# that a transcript-space NH tag would conflate it with.
mcols(mapped)$NH <- nh[mapped$xHits]
message("Genomic NH distribution among projected rows: ",
        paste(capture.output(print(table(cut(mcols(mapped)$NH, c(0,1,2,4,8,Inf),
                                              labels = c("1","2","3-4","5-8",">8"))))),
              collapse = " | "))
names(mapped) <- read_names[mapped$xHits]
mapped$xHits <- NULL
mapped$transcriptsHits <- NULL
# mapToTranscripts() returns transcript-space coordinates that already run
# 5' -> 3' along the transcript, but sets strand to the *genomic* strand of
# the transcript (~half "-" for genes on the reverse strand). Ribostan's
# get_psite_gr() then does resize(1, "start") before applying the P-site
# offset: for strand="-" that jumps to the numerical end of the range
# (Bioconductor's biological-5' convention), which is the wrong end in
# already-oriented transcript space and destroys periodicity for every
# reverse-strand gene. The damage is frame-silent when (readlen-1) %% 3 == 0
# (e.g. 28 nt stays fine) and catastrophic otherwise (26/27/29/30 nt drop
# to ~chance) - exactly the per-length pattern seen in the HEK293 genome-
# route diagnostics. Transcriptome BAMs never hit this: read_ribobam()
# keeps only isMinusStrand=FALSE, and STAR transcript alignments are all
# "+" relative to the transcript FASTA. Force "+" here so transcript-space
# start() remains the 5' end under resize().
strand(mapped) <- "+"
mapped <- sort(mapped)
# Keep the true read names (e.g. "SRR2433794.28659155") rather than
# collapsing them to a compact integer code via match(names, unique(names)):
# the periodicity pipeline only needs name equality to be preserved (which a
# match()-based code also satisfies), but downstream forensics that
# cross-reference these reads against the transcriptome BAM's own read names
# (e.g. "is this genome-projected read also alignable in transcript space at
# all") need the real identifiers.

saveRDS(mapped, out_rds)
message("Saved ", length(mapped), " transcript-space reads to ", out_rds)
