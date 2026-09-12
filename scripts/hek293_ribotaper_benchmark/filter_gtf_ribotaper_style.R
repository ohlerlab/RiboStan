#!/usr/bin/env Rscript
# Build a Ribostan-loadable GTF using RiboTaper's default GENCODE v19 filters
# (create_annotations_files.bash: use_ccdsid=true, use_appris=true).
args <- commandArgs(trailingOnly = FALSE)
script_path <- sub("--file=", "", args[grep("--file=", args)])
bench_dir <- if (length(script_path)) dirname(script_path) else "."
workdir <- {
  tra <- commandArgs(trailingOnly = TRUE)
  if (length(tra) >= 1) tra[[1]] else "."
}
setwd(workdir)

source(file.path(bench_dir, "hek293_gtf_config.R"), local = TRUE)

rebuild <- as.logical(Sys.getenv("HEK293_REBUILD_GTF", "FALSE"))
cfg <- resolve_hek293_gtf(workdir, rebuild = rebuild)
if (cfg$built) {
  message("Wrote ", cfg$path)
  print(cfg$summary)
} else if (isTRUE(cfg$explicit)) {
  message("Using explicit HEK293_GTF_NAME=", cfg$name, " (no filter step)")
} else {
  message("Reusing existing filtered GTF: ", cfg$name)
  if (file.exists(cfg$summary_path)) {
    print(read.delim(cfg$summary_path, stringsAsFactors = FALSE))
  }
}
