# HEK293 RiboTaper benchmark (Calviello et al. 2016)

Reference dataset where RiboTaper **did** call translated uORFs.

## Source

| Item | Accession |
|------|-----------|
| Study | [GSE73136](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE73136) |
| Ribo-seq | [SRR2433794](https://www.ncbi.nlm.nih.gov/sra/SRR2433794) (HEK293, ~32M reads, 27–30 nt footprints) |
| Published ORFs | `GSM1887643_HEK293_filtered_ORFs.bed.gz` (GEO supplementary) |
| Paper | Calviello et al., *Nature Methods* 2016 ([nmeth.3688](https://doi.org/10.1038/nmeth.3688)) |
| Annotation | Gencode v19 / GRCh37 (same as original analysis) |

## RiboTaper uORF calls (ground truth)

The filtered BED contains **908 uORF intervals** across **701 transcripts**
(category `uORF` in the BED name field, e.g.
`ENST00000428771.2_9_66;uORF;protein_coding`).

These are translation-supported: periodic RPF signal + ORF reconstruction,
not sequence-only candidates.

## Ribostan comparison workflow (RiboTaper-matched — use this)

The fair comparison uses **RiboTaper's own read path**, not Ribostan's
transcriptome BAM default:

1. `download_data.sh` — Gencode v19, GEO BED, SRR2433794 FASTQ
2. `filter_rrna_bowtie1.sh` — Bowtie1 `--best` rRNA pre-filter (GSM1887643 notes)
3. `align_genome_ribotaper_style_rrnafilt.sh` — genome STAR + primary-only
4. `convert_genome_to_txspace.R` — project to transcript space (`strand="+"` fix;
   keep spliced CIGARs; dedupe `anno$exonsgrl`)
5. `compare_hek293_prep.R` / `compare_hek293_report.R` — defaults now match
   RiboTaper (CCDS/APPRIS GTF, `metagene_peak` offsets, `ribotaper` nested rule,
   raw p<0.05, in-frame gate 0.5, `genome_primary` reads)
6. `diagnose_differences.R` — per-uORF miss/excess classification

One command (skips steps whose outputs already exist):

```bash
bash scripts/hek293_ribotaper_benchmark/run_ribotaper_matched.sh /path/to/workdir
```

On EC2:

```bash
bash scripts/ec2_hek293_comparison_conda.sh
```

Legacy transcriptome-only route (`align_transcriptome.sh` +
`HEK293_READ_SOURCE=transcriptome_bam`) remains for Ribostan-native analysis
but is **not** the RiboTaper agreement benchmark.

## Exact parity baseline (run before comparisons)

Before changing offsets, reads, annotations, or filtering, validate the two
periodicity paths on literally shared inputs:

```bash
Rscript scripts/hek293_ribotaper_benchmark/validate_exact_parity.R /path/to/workdir
```

This imports one RiboTaper P-site track, uses the published RiboTaper uORFs,
constructs each ORF coverage vector once, and passes those vectors through the
original RiboTaper calculation and RiboStan's `ftest_orfs` path. The script
stops on any p-value, frame-fraction, or call disagreement. A successful run
prints `EXACT PARITY`.

Only after this reaches 100% should experiments replace one component at a
time, in this order: candidate ORFs, P-site construction with fixed read
lengths/offsets, offset estimation, then other RiboStan-native choices.

## Methodology-matching knobs

`compare_hek293_prep.R` accepts environment variables to align Ribostan's
periodicity filtering with RiboTaper's own published methodology
(Calviello et al. 2016), so the two can be compared on identical
annotation/data with as few incidental differences as possible:

| Variable | Values | RiboTaper-matching setting |
|---|---|---|
| `HEK293_GTF_SOURCE` | GTF filename in workdir | **`gencode.v19.annotation.gtf`** (full GENCODE v19 — what RiboTaper passes to `create_annotations_files.bash`) |
| `HEK293_USE_CCDS` | `TRUE`/`FALSE` | **`TRUE`** (`use_ccdsid=true`: CCDS genes = genes with a `ccdsid`-tagged CDS row) |
| `HEK293_USE_APPRIS` | `TRUE`/`FALSE` | **`TRUE`** (`use_appris=true`: one isoform per CCDS gene — `appris_principal` first, else any APPRIS tag, else all isoforms) |
| `HEK293_GTF_NAME` | GTF filename in workdir | *(unset by default)* — auto-builds `gencode.v19.ccds_appris.gtf` from source + flags above; set to skip auto-build |
| `HEK293_REBUILD_GTF` | `TRUE`/`FALSE` | `FALSE` — force regeneration of the filtered GTF |
| `HEK293_FDR_CORRECT` | `TRUE`/`FALSE` | `FALSE` (RiboTaper thresholds on the raw per-ORF p-value, not a BH q-value) |
| `HEK293_NESTED_RULE` | `density`/`ribotaper` | `ribotaper` (RiboTaper's own upstream-ATG-with-signal / "max_P_sites" rule for resolving ATGs sharing a stop codon - see `select_ribotaper_nested_uorfs`) |
| `HEK293_MULTIMAP` | `all`/`primary`/`unique` | `all` (best empirical match; `primary`/`unique` lose real periodicity signal without a precision benefit) |
| `HEK293_MIN_FRAME_FRAC` | `0`-`1` | `0.5` (RiboTaper's `Phase_P_sites_frame/P_sites_sum > 0.5` gate) |
| `HEK293_OFFSET_METHOD` | `cds_incl`/`metagene_peak` | `metagene_peak` (RiboTaper's own visual metagene-peak P-site offset calibration, reproduced programmatically) |
| `HEK293_APPLY_PHASE_SHIFT` | `TRUE`/`FALSE` | **`FALSE`** (RiboTaper uses metagene offset only; Ribostan's default CDS phase-shiftdf reranking is off for this benchmark) |
| `HEK293_MULTIMAP_QC` | `TRUE`/`FALSE` | see caveat below |
| `HEK293_N_CORES` | integer | parallelizes the (slow) multitaper periodicity test only |
| `HEK293_READ_SOURCE` | `transcriptome_bam`/`genome_primary` | **`genome_primary`** (default; RiboTaper's genome STAR + primary-only path) |
| `HEK293_PROJECTED_RDS` | RDS filename in workdir | **`rpfs_genome_projected_rrnafilt_matched.rds`** (rRNA-filtered genome BAM projected with strand fix + full CIGAR) |

Caveat on `HEK293_MULTIMAP_QC`: this reproduces RiboTaper's post-hoc
`flag_multimapping_uorfs` coverage QC (drop ORFs whose read coverage is
dominated by multi-mapping reads and/or too sparse overall), but RiboTaper
calibrated its 0.3 ratio threshold on **genomic** multi-mapping (ambiguous
between distinct genomic loci - repeats/paralogs). Ribostan aligns to the
**transcriptome**, where a read landing on several isoforms of the *same*
gene's shared 5' UTR also counts as multi-mapping (`NH > 1`), even though
it is unambiguous at the genomic level. Applying RiboTaper's threshold
directly here flags the large majority of real uORFs (many have
`pct_covered_unique == 0` simply because every isoform shares the
sequence), so treat this knob as diagnostic rather than a recommended
default until multi-mapping is redefined at the genomic-locus level.

## Genome-alignment experiment (`HEK293_READ_SOURCE=genome_primary`)

RiboTaper's own published pipeline (per GSM1887643's "Data processing" notes)
aligns to the **genome** (STAR + Gencode v19 splice junctions), then keeps
only the primary alignment per read (`samtools view -F 0x100`) - i.e. exactly
one genomic locus per read, discarding true multi-mappers beyond
`--outFilterMultimapNmax 8` entirely. `align_genome_ribotaper_style.sh` +
`convert_genome_to_txspace.R` reproduce this: build a genome STAR index,
align with RiboTaper's filters (`--alignEndsType EndToEnd
--outFilterMismatchNmax 4 --outFilterMultimapNmax 8
--outFilterIntronMotifs RemoveNoncanonicalUnannotated --alignSJoverhangMin
500`), keep primary-only, then project each read's genomic span into every
transcript whose exons it falls in (`GenomicFeatures::mapToTranscripts`
against `anno$exonsgrl`) - the transcript-space analogue of `multimap="all"`.
Setting `HEK293_READ_SOURCE=genome_primary` makes `compare_hek293_prep.R`
load this pre-projected read set (`rpfs_genome_projected.rds`) instead of
calling `get_readgr()` on the transcriptome BAM.

**Result (before the strand fix below): this made recall worse, not
better** (RiboTaper uORFs recovered: 27.4% vs 58.0% for the
transcriptome/`multimap="all"` route, both otherwise using identical
settings). After forcing `strand(mapped) <- "+"` post-projection (see root
cause below), genome-route recall rises to **48.7%**. The pre-fix
genome-alignment route pulled in ~3.8x more P-sites per candidate on
average, but the extra volume was disproportionately non-periodic: for
candidates that were significant via the transcriptome route but lost
significance here, median `pct_inframe` dropped from 0.68 to 0.30 - right
at the 1-in-3 random-chance floor - while median depth *tripled*.
Two specific causes were ruled out directly:
- **Isoform-redundant candidate annotation** (the same biological uORF
  double-counted under >1 `transcript_id`): only 80/27,998 candidate loci
  are duplicated this way; excluding them changes nothing.
- **True genomic multi-mapping/repeat noise**: restricting to reads that
  were already unique at the genomic level (`NH==1`, 88.8% of the projected
  signal) barely moves the needle (26.7% vs 27.4% recall) - the dilution is
  present even among uniquely-mapping reads, so it isn't a multi-mapping
  problem Ribostan's existing `flag_multimapping_uorfs` QC could fix.

**rRNA pre-filtering was tested directly and ruled out.** GSM1887643's GEO
"Data processing" notes (verified by re-fetching the record, not assumed)
state the real pipeline runs Bowtie1 (`--best`) against rRNA first and keeps
only the unaligned reads before genome alignment - a step this experiment
had omitted. `filter_rrna_bowtie1.sh` reproduces it (human 45S/5S nuclear +
mitochondrial 12S/16S rRNA reference, Bowtie1 v1.3.1 `--best`), and
`align_genome_ribotaper_style_rrnafilt.sh` +
`convert_genome_to_txspace.R` (`HEK293_PROJECTED_RDS=rpfs_genome_projected_rrnafilt.rds`)
re-run the rest of the genome route on the filtered reads
(`HEK293_OUT_SUFFIX=rrnafilt`). Of 31.97M raw reads, 8.3% (2.65M) aligned to
rRNA and were removed. **Recall was unchanged: 27.6% vs 27.4% without the
filter** (251/908 vs 252/908 RiboTaper uORFs recovered). The significant-call
sets from the two genome routes overlap almost perfectly (527/534 candidates
agree; only 7 flip each way), and the lost-vs-kept pattern from before is
identical: for candidates significant via the transcriptome route, both the
"kept" and "lost" (in the genome route) groups see ~3x more P-sites, but only
the "lost" group's `pct_inframe` collapses (0.72 -> 0.28, chance level) while
the "kept" group's does not (0.73 -> 0.73). So the extra non-periodic depth
recruited by genome alignment is not rRNA - rRNA reads are too rare (8.3%)
and too evenly distributed between the kept/lost sets to explain a
recall-halving effect concentrated in specific 5' UTR candidates.

**A real bug was found and fixed in the diagnostic script - but it turned out
not to be the explanation.** `anno$exonsgrl` (from `load_annotation()`) is
*not* transcript-unique: it has one entry per **ORF** (one CDS + every uORF
candidate on that transcript), not one per transcript, because
`load_annotation()` indexes a transcript-keyed exon list with a per-ORF
`transcript_id` vector for internal use by `get_trspace_cds()`, then exposes
that replicated object as the public `anno$exonsgrl` field (its length is
pinned at 3175 for the chr22 test fixture by
`test-annotation.R`, so this is long-standing, tested package behavior, not
something to change unilaterally). `convert_genome_to_txspace.R` used
`anno$exonsgrl` as a whole-list `mapToTranscripts()` projection target rather
than doing by-name lookups, so every read landing on a transcript with *N*
annotated ORFs was broadcast into *N* near-identical output rows. Confirmed
directly: HEK293 candidate `ENST00000598463.1_6` (on a transcript with 1 CDS
+ 9 uORF candidates = 10 ORFs) had exactly 10x as many P-sites as unique
reads. Globally, 27,998 of `anno$exonsgrl`'s 46,611 entries were excess
duplicates - exactly equal to the total uORF candidate count, confirming the
"one extra copy per uORF" mechanism.

The fix (deduplicate `anno$exonsgrl` by transcript name before calling
`mapToTranscripts()`, done in `convert_genome_to_txspace.R`) cut total
projected P-sites roughly in half (33.7M -> 15.6M) and brought the
read-count-to-P-site ratio back to a sane ~1.006x (matching genuine isoform
sharing, not artificial broadcast). **But recall was statistically
unchanged: 27.4% both before and after the fix** (249/908 vs 251/908
RiboTaper uORFs recovered), and re-running the corrected per-read-length CDS
frame-purity check on the fixed P-sites reproduced the same degradation
pattern as before the fix (25nt 47%, 26nt 40%, 27nt 42%, 28nt 72%, 29nt 46%,
30nt 48% frame-0, vs the transcriptome route's 47/65/58/72/85/72% for the
same lengths). This **definitively rules out** the `exonsgrl` duplication bug
as the cause of the genome-route's poor recall: it was real and worth fixing,
but the underlying non-periodic dilution is present in the deduplicated data
too, at essentially the same magnitude.

Further hypotheses tested after the `exonsgrl` fix (all on the
`rrnafilt_fixed_names` projected set):

- **Genome-only "spurious locus" rescue**: 99.7% of genome-route CDS read
  names are also present somewhere in the transcriptome BAM. Ruled out.
- **Wrong-gene primary-locus assignment**: 99.7% of those reads are
  gene-consistent with their transcriptome-BAM alignments. Ruled out.
- **Per-nucleotide placement disagreement** (genome+`mapToTranscripts` vs
  transcriptome STAR, same `(read_name, transcript_id)`): **100.00% exact
  position agreement** after de-fanning join duplicates (15,183,829 / 
  15,184,060 matched pairs; only 231 residual disagreements, mostly ±1-2
  nt). Ruled out as the cause of the periodicity gap.

**Root cause found: `mapToTranscripts()` strand vs transcript-space
coordinates.** Projected ranges already use 5'→3' transcript coordinates,
but `mapToTranscripts()` sets `strand` to the transcript's *genomic* strand
(~50% `"-"` for reverse-strand genes). Ribostan's `get_psite_gr()` then does
`resize(1, "start")` before applying the P-site offset. For `strand="-"`
that jumps to the numerical *end* of the range (Bioconductor's
biological-5' convention) - the wrong end in already-oriented transcript
space - and destroys periodicity for every reverse-strand gene.

This predicts (and matches) the observed per-length pattern exactly:
`(readlen - 1) %% 3 == 0` lengths are frame-silent under a
start↔end flip (28 nt stayed at ~72% frame-0, matching the transcriptome
route), while 26/27/29/30 nt collapsed toward chance. Raw `start()`-based
frames before `resize` are ~75% frame-0 on *both* strands (matching the
transcriptome route); after `resize(1, "start")`, `"+"` stays at 75% and
`"-"` drops to ~23%. Transcriptome BAMs never hit this: `read_ribobam()`
keeps only `isMinusStrand=FALSE`, and STAR transcript alignments are all
`"+"` relative to the transcript FASTA.

**Fix** (in `convert_genome_to_txspace.R`): force `strand(mapped) <- "+"`
after projection so transcript-space `start()` remains the 5' end under
`resize()`. Re-run with
`HEK293_PROJECTED_RDS=rpfs_genome_projected_rrnafilt_strandfix.rds` /
`HEK293_OUT_SUFFIX=rrnafilt_strandfix` recovered **48.7%** of RiboTaper
uORFs (442/908; RiboTaper-matched definition) vs **27.4%** before the
strand fix and **58.0%** for the transcriptome/`multimap="all"` route.
That closes most of the genome-route gap; the residual ~9 pp vs
transcriptome likely reflects remaining retention differences (genome
`MultimapNmax=8` + primary-only vs transcriptome `MultimapNmax=20` +
`multimap="all"`, and `simpleCigar=TRUE` dropping ~13.6% junction-spanning
genome alignments that the transcriptome route keeps as contiguous
transcript hits) - not yet quantified separately. Prefer the default
transcriptome/`multimap="all"` route for real analysis; the genome route
remains a diagnostic knob, now with a known, fixed strand footgun.

**Note for anyone else consuming `anno$exonsgrl` directly:** if you use it
as a whole-list target for `mapToTranscripts()`, `extractTranscriptSeqs()`,
or similar functions that iterate every list element (rather than doing
`anno$exonsgrl[transcript_id]`-style by-name lookups, which silently return
just the first match and are unaffected), deduplicate first:
`anno$exonsgrl[!duplicated(names(anno$exonsgrl))]`.

## Notes

- RNA-seq for RiboTaper normalization was GSM1306496 (external study); we use **published uORF BED** as reference rather than re-running RiboTaper.
- Expect several hours on EC2 (download ~4 GB annotation + ~1 GB FASTQ, STAR alignment, full-genome uORF scan).
