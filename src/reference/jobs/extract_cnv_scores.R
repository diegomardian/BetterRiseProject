# Recover per-cell CNV scores from a finished inferCNV run. W1, weeks 2-3.
#
#   Rscript src/reference/jobs/extract_cnv_scores.R <out_dir> [<out_dir> ...]
#
# For runs that completed BEFORE run_infercnv() started writing cnv_scores.csv
# itself. inferCNV's no_plot = TRUE skips the step that writes
# infercnv.observations.txt, so those runs finished with the answer inside
# run.final.infercnv_obj and nothing readable beside it.
#
# Recovering is minutes; re-running the inference is hours per patient. There is
# no reason to pay the second price for the first mistake.
#
# The score is the mean squared deviation from 1 across genes — the same
# quantity run_infercnv() now computes inline, and what call_malignancy()
# thresholds. Squared rather than absolute: a cell with a few large deviations
# is more plausibly aneuploid than one with many tiny ones.

suppressPackageStartupMessages(library(infercnv))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("usage: Rscript extract_cnv_scores.R <out_dir> [<out_dir> ...]")
}

for (dir in args) {
  final <- file.path(dir, "run.final.infercnv_obj")
  if (!file.exists(final)) {
    cat("SKIP", dir, "- no run.final.infercnv_obj (did the run finish?)\n")
    next
  }
  out <- file.path(dir, "cnv_scores.csv")
  if (file.exists(out)) {
    cat("SKIP", dir, "- cnv_scores.csv already present\n")
    next
  }

  obj <- readRDS(final)
  expr <- obj@expr.data
  scores <- colMeans((expr - 1)^2)
  write.csv(
    data.frame(cell = names(scores), cnv_score = as.numeric(scores)),
    file = out, row.names = FALSE
  )
  cat("wrote", out, "-", length(scores), "cells,", nrow(expr), "genes\n")
}
