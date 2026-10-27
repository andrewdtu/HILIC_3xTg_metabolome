#!/usr/bin/env bash
# Reproduce Figure 6 from the supplied protein-normalized table. Run from the repository root.
set -euo pipefail
Rscript R/01_factorial_limma.R none   # protein-normalized (as supplied)
Rscript R/01_factorial_limma.R pqn    # protein + global PQN (primary)
# Rscript R/02_pathway_sets.R         # optional: rebuild annotation/pathway_sets_all_dbs.rds (needs annotation/db/, see README)
Rscript R/03_enrichment.R pqn
Rscript R/04_figure6.R
Rscript -e "suppressPackageStartupMessages({library(limma); library(dplyr); library(tidyr); library(ggplot2); library(patchwork); library(ragg)}); writeLines(capture.output(sessionInfo()), \"results/sessionInfo.txt\")"
