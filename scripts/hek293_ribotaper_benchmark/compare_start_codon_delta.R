#!/usr/bin/env Rscript
# Compare AUG-only calls with the isolated addition of CTG starts.
args <- commandArgs(trailingOnly = TRUE)
workdir <- if (length(args)) args[[1L]] else "."
setwd(workdir)

suppressPackageStartupMessages({
  library(Ribostan)
  library(GenomicRanges)
  library(rtracklayer)
})

record_path <- Sys.getenv(
  "HEK293_START_CODON_RECORDS",
  "psite_records_general_26_30.tsv"
)
output_path <- Sys.getenv(
  "HEK293_START_CODON_OUTPUT",
  "ribostan_start_codon_delta.tsv"
)
n_cores <- as.integer(Sys.getenv("HEK293_N_CORES", "1"))
reference_path <- Sys.getenv(
  "HEK293_PARITY_UORFS",
  "ribotaper_hek293_uorfs.tsv"
)
gtf_path <- Sys.getenv("HEK293_FULL_GTF", "gencode.v19.annotation.gtf")
fasta_path <- Sys.getenv("HEK293_FASTA", "GRCh37.p13.genome.fa")

reference <- read.delim(reference_path, stringsAsFactors = FALSE)
target_transcripts <- unique(reference$transcript_id)
gtf <- import(gtf_path, format = "gtf")
target_genes <- unique(as.character(
  gtf$gene_id[as.character(gtf$transcript_id) %in% target_transcripts]
))
keep <- as.character(gtf$transcript_id) %in% target_transcripts |
  (is.na(gtf$transcript_id) & as.character(gtf$gene_id) %in% target_genes)
target_gtf <- tempfile(fileext = ".gtf")
on.exit(unlink(target_gtf), add = TRUE)
export(gtf[keep], target_gtf, format = "gtf")
rm(gtf)

records <- read.delim(record_path, stringsAsFactors = FALSE)
records <- records[
  !is.na(records$position_ribostan) &
    records$read_length %in% 26:30,
  ,
  drop = FALSE
]
records <- records[!duplicated(records[
  c("read_name", "transcript_id", "position_ribostan")
]), , drop = FALSE]
raw_psites <- GRanges(
  records$transcript_id,
  IRanges(records$position_ribostan, width = 1L),
  strand = "+"
)
names(raw_psites) <- records$read_name

evaluate_scenario <- function(scenario, start_codons, candidate_path,
                              nested_rule = "ribotaper",
                              input_psites = raw_psites) {
  anno <- load_annotation(
    target_gtf,
    fasta_path,
    add_uorfs = TRUE,
    findUORFs_args = list(startCodon = start_codons)
  )
  candidates <- read.delim(candidate_path, stringsAsFactors = FALSE)
  uorf_ids <- names(anno$uORF)[anno$uORF]
  orfs <- anno$trspacecds[uorf_ids]
  overlaps <- findOverlaps(input_psites, orfs, ignore.strand = TRUE)
  psites <- input_psites[queryHits(overlaps)]
  mcols(psites)$orf <- names(orfs)[subjectHits(overlaps)]

  tested <- periodicity_filter_uORFs(
    psites,
    anno,
    remove = FALSE,
    n_cores = n_cores,
    nested_rule = nested_rule,
    fdr_correct = TRUE,
    min_frame_frac = 0.5,
    min_psites = 6L
  )
  metadata <- as.data.frame(mcols(tested$trspacecds[uorf_ids]))
  result <- candidates[match(uorf_ids, candidates$orf_id), , drop = FALSE]
  result$scenario <- scenario
  for (field in c(
    "p.value", "q.value", "n_psites", "pct_inframe", "orf_locus_id",
    "locus_representative"
  )) {
    result[[field]] <- metadata[[field]]
  }
  supported <- !is.na(result$n_psites) &
    result$n_psites >= 6L &
    !is.na(result$pct_inframe) &
    result$pct_inframe > 0.5
  result$call_raw <- supported &
    !is.na(result$p.value) &
    result$p.value < 0.05
  result$call_bh <- supported &
    !is.na(result$q.value) &
    result$q.value < 0.05
  result
}

scenario_results <- list(
  evaluate_scenario(
    "ATG",
    "ATG",
    "ribostan_candidates_atg.tsv"
  ),
  evaluate_scenario(
    "ATG_CTG",
    c("ATG", "CTG"),
    "ribostan_candidates_atg_ctg.tsv"
  ),
  evaluate_scenario(
    "ATG_DENSITY",
    "ATG",
    "ribostan_candidates_atg.tsv",
    nested_rule = "density"
  )
)
if (identical(
  toupper(Sys.getenv("HEK293_INCLUDE_UNIQUE_SCENARIO", "FALSE")),
  "TRUE"
)) {
  bam <- Sys.getenv(
    "HEK293_GENOME_BAM",
    "hek293_genome_rrnafilt_primary.bam"
  )
  bam_records <- Rsamtools::scanBam(
    bam,
    param = Rsamtools::ScanBamParam(what = "qname", tag = "NH")
  )[[1L]]
  nh_by_name <- as.integer(bam_records$tag$NH)
  names(nh_by_name) <- bam_records$qname
  unique_psites <- raw_psites[nh_by_name[names(raw_psites)] == 1L]
  scenario_results[[length(scenario_results) + 1L]] <- evaluate_scenario(
    "ATG_UNIQUE",
    "ATG",
    "ribostan_candidates_atg.tsv",
    input_psites = unique_psites
  )
}
results <- do.call(rbind, scenario_results)

write.table(
  results,
  output_path,
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)
loci <- results[
  !is.na(results$locus_representative) & results$locus_representative,
  ,
  drop = FALSE
]
print(with(loci, table(scenario, start_codon)))
print(with(loci, table(scenario, start_codon, call_raw)))
print(with(loci, table(scenario, start_codon, call_bh)))
message("Wrote ", output_path)
