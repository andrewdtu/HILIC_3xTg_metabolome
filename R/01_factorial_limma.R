# 2x2x2 (Genotype x Sex x Age) factorial limma analysis of the 3xTg hippocampal metabolome.
# Input is already normalized to BCA protein (ng / ug protein). Usage: Rscript R/01_factorial_limma.R {none|pqn}
suppressPackageStartupMessages({library(limma); library(dplyr); library(tidyr)})
args <- commandArgs(TRUE); NORM <- if (length(args)) args[1] else "none"  # "none" (BCA protein, as supplied) or "pqn" (protein + PQN)
stopifnot(NORM %in% c("none", "pqn"))
cat("Normalization:", NORM, "\n")
d <- read.csv("data/3xTg_hippocampus_metabolites_protein_normalized.csv", check.names = FALSE)
names(d)[1] <- "Sample"
meta <- d[, 1:2]
meta$Genotype <- factor(ifelse(grepl("3xTg", meta$Sample), "Tg", "WT"), levels = c("WT", "Tg"))
meta$Sex <- factor(ifelse(grepl("-F-", meta$Sample), "F", "M"))
meta$Age <- factor(ifelse(grepl("^12-", meta$Sample), "12mo", "3mo"), levels = c("3mo", "12mo"))
meta$cell <- factor(paste(meta$Genotype, meta$Sex, meta$Age, sep = "_"))
print(table(meta$cell))

X <- as.matrix(d[, -(1:2)])
# zeros -> half of the per-metabolite minimum positive value, then log2
X <- apply(X, 2, function(v) { v[v <= 0] <- min(v[v > 0]) / 2; v })
pqn <- function(X) { ref <- apply(X, 2, median); apply(sweep(X, 2, ref, "/"), 1, median) }  # probabilistic quotient normalization
if (NORM == "pqn") {
  q <- pqn(X); X <- X / q; meta$pqn_factor <- q
}
Y <- t(log2(X)); colnames(Y) <- meta$Sample

design <- model.matrix(~ 0 + cell, meta); colnames(design) <- levels(meta$cell)
c_ <- function(g, s, a) paste(g, s, a, sep = "_")
G <- function(s, a) paste0("(", c_("Tg", s, a), " - ", c_("WT", s, a), ")")
cons <- makeContrasts(
  Geno_F3   = paste(G("F","3mo")),
  Geno_F12  = paste(G("F","12mo")),
  Geno_M3   = paste(G("M","3mo")),
  Geno_M12  = paste(G("M","12mo")),
  GxA_F     = paste(G("F","12mo"), "-", G("F","3mo")),
  GxA_M     = paste(G("M","12mo"), "-", G("M","3mo")),
  GxA_avg   = paste("(", G("F","12mo"), "+", G("M","12mo"), "-", G("F","3mo"), "-", G("M","3mo"), ")/2"),
  GxS_3     = paste(G("F","3mo"), "-", G("M","3mo")),
  GxS_12    = paste(G("F","12mo"), "-", G("M","12mo")),
  GxSxA     = paste("(", G("F","12mo"), "-", G("F","3mo"), ") - (", G("M","12mo"), "-", G("M","3mo"), ")"),
  Geno_main = paste("(", G("F","3mo"), "+", G("F","12mo"), "+", G("M","3mo"), "+", G("M","12mo"), ")/4"),
  levels = design)
colnames(cons) <- c("Geno_F3","Geno_F12","Geno_M3","Geno_M12","GxA_F","GxA_M","GxA_avg","GxS_3","GxS_12","GxSxA","Geno_main")
fit <- eBayes(contrasts.fit(lmFit(Y, design), cons))

res <- bind_rows(lapply(colnames(cons), function(k)
  topTable(fit, coef = k, number = Inf, sort.by = "none") %>%
    tibble::rownames_to_column("Metabolite") %>% mutate(contrast = k)))
OUT <- file.path("results", NORM); dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
write.csv(res, file.path(OUT, "limma_all_contrasts_long.csv"), row.names = FALSE)

wide <- res %>% select(Metabolite, contrast, logFC, P.Value, adj.P.Val) %>%
  pivot_wider(names_from = contrast, values_from = c(logFC, P.Value, adj.P.Val), names_glue = "{contrast}_{.value}")
write.csv(wide, file.path(OUT, "limma_all_contrasts_wide.csv"), row.names = FALSE)

# omnibus F-tests for factorial terms (sum-to-zero coding), per metabolite
options(contrasts = c("contr.sum", "contr.poly"))
dsg2 <- model.matrix(~ Genotype * Sex * Age, meta)
fit2 <- eBayes(lmFit(Y, dsg2))
terms <- colnames(dsg2)[-1]
omni <- sapply(terms, function(t) fit2$p.value[, t])
omni_fdr <- apply(omni, 2, p.adjust, "BH")
colnames(omni) <- paste0("p_", colnames(omni)); colnames(omni_fdr) <- paste0("fdr_", colnames(omni_fdr))
write.csv(data.frame(Metabolite = rownames(Y), omni, omni_fdr, check.names = FALSE), file.path(OUT, "factorial_terms_pvals.csv"), row.names = FALSE)

cat("\n# metabolites p<0.05 / FDR<0.1 per contrast\n")
print(res %>% group_by(contrast) %>% summarise(p05 = sum(P.Value < .05), up = sum(P.Value < .05 & logFC > 0),
        down = sum(P.Value < .05 & logFC < 0), fdr10 = sum(adj.P.Val < .1), fdr05 = sum(adj.P.Val < .05)))
cat("\n# factorial terms: n p<0.05, n FDR<0.1\n")
print(rbind(p05 = colSums(omni < .05), fdr10 = colSums(omni_fdr < .1)))
saveRDS(list(Y = Y, meta = meta, fit = fit, res = res, wide = wide), file.path(OUT, "fit.rds"))
