# Enrichment across KEGG, Reactome, SMPDB and RefMet chemical classes. Usage: Rscript R/03_enrichment.R [pqn]
# Directional ORA (hypergeometric on p < 0.05 hits, up and down separately; used in Figure 6C) and camera
# (competitive, correlation-adjusted). Background = measured metabolites annotated in each resource; sets with >= 3 members.
suppressPackageStartupMessages({library(limma); library(dplyr); library(tidyr)})
sets <- readRDS("annotation/pathway_sets_all_dbs.rds"); sets$Reactome_all_levels <- NULL
rm_ <- read.csv("annotation/refmet_classes.csv")
sets$RefMet_main_class <- rm_ %>% transmute(id = main_class, name = main_class, Metabolite)
sets$RefMet_sub_class <- rm_ %>% transmute(id = sub_class, name = sub_class, Metabolite)

NORM <- if (length(commandArgs(TRUE))) commandArgs(TRUE)[1] else "pqn"
r <- readRDS(sprintf("results/%s/fit.rds", NORM)); Y <- r$Y; design <- r$fit$design; cons <- r$fit$contrasts
out <- list()
for (d in names(sets)) {
  S <- sets[[d]] %>% filter(Metabolite %in% rownames(Y)) %>% group_by(id, name) %>% summarise(members = list(sort(unique(Metabolite))), n = length(unique(Metabolite)), .groups = "drop") %>%
    filter(n >= 3) %>% distinct(members, .keep_all = TRUE)  # drop sets with identical measured members (sorted, so order does not matter; common in SMPDB)
  bg <- unique(sets[[d]]$Metabolite); Yd <- Y[bg, ]  # background = measured metabolites annotated in this resource
  idx <- setNames(lapply(S$members, function(m) match(m, rownames(Yd))), S$id)
  for (k in colnames(cons)) {
    cam <- camera(Yd, idx, design, contrast = cons[, k], inter.gene.cor = NA) %>% tibble::rownames_to_column("id")
    tt <- r$res %>% filter(contrast == k, Metabolite %in% bg)
    ora <- function(dir) { hits <- tt$Metabolite[tt$P.Value < .05 & sign(tt$logFC) == dir]
      sapply(S$members, function(m) { x <- sum(m %in% hits); if (x == 0) 1 else phyper(x - 1, length(m), nrow(tt) - length(m), length(hits), lower.tail = FALSE) }) }
    out[[paste(d, k)]] <- tibble(database = d, contrast = k, id = S$id, name = S$name, n = S$n,
      mean_logFC = sapply(S$members, function(m) mean(tt$logFC[tt$Metabolite %in% m])),
      hits = sapply(S$members, function(m) { h <- tt %>% filter(Metabolite %in% m, P.Value < .05); paste0(h$Metabolite, ifelse(h$logFC > 0, "↑", "↓"), collapse = ", ") }),
      ora_up_p = ora(1), ora_down_p = ora(-1)) %>%
      mutate(ora_up_fdr = p.adjust(ora_up_p, "BH"), ora_down_fdr = p.adjust(ora_down_p, "BH")) %>%
      left_join(cam %>% select(id, camera_dir = Direction, camera_p = PValue, camera_fdr = FDR), by = "id")
  }
}
res <- bind_rows(out) %>% arrange(database, contrast, camera_p)
write.csv(res, sprintf("results/enrichment_multi_db_%s.csv", NORM), row.names = FALSE)
cat("Sets tested per database:\n"); print(res %>% distinct(database, id) %>% count(database))
cat("\nFDR < 0.1 (camera or ORA):\n"); options(width = 250)
print(as.data.frame(res %>% filter(camera_fdr < .1 | ora_up_fdr < .1 | ora_down_fdr < .1) %>%
  transmute(database, contrast, name = substr(name, 1, 45), n, mFC = round(mean_logFC, 2), cam = camera_dir, cam_fdr = signif(camera_fdr, 2),
            oraU = signif(ora_up_fdr, 2), oraD = signif(ora_down_fdr, 2), hits = substr(hits, 1, 80))), right = FALSE)
