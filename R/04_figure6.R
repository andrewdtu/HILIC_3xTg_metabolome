# Figure 6: 3xTg hippocampal metabolome (protein + global PQN unless noted). Usage: Rscript R/04_figure6.R
#  A (left half): log2FC heatmap, all 92 metabolites, protein-only vs protein + PQN, genotype x age
#  Right half, females only:
#  B  number of metabolites up / down in 3xTg (pale: nominal p < 0.05, solid: FDR < 0.05)
#  C  directional ORA pathway enrichment (p = best of up/down; rows: sets with p < 0.05 in any of the 7 comparisons)
#  D  lysine / 2-HG metabolites by age (mean +/- SE, mice as dots), log10(abundance + 1)
#  E  lysine degradation pathway (after Tan, Chrysopoulou & Rinschen 2023), F 3 mo and F 12 mo
# * FDR < 0.05, # nominal p < 0.05 throughout.
suppressPackageStartupMessages({library(dplyr); library(tidyr); library(ggplot2); library(patchwork)})
source("R/theme_pew_serif.R")  # lab house theme (Libre Franklin / Source Serif 4 fonts)
NORM <- "pqn"; CMP <- "none"  # panel A compares protein + PQN with protein-only
nlab <- c(none = "Protein-normalized", pqn = "Protein + PQN")
f <- readRDS(sprintf("results/%s/fit.rds", NORM))
up <- "#e34948"; dn <- "#2a78d6"; mid <- "#f0efec"
BS <- 15                                    # base text size (pt) for every panel
gs <- BS / ggplot2::.pt                     # geom text size matching base
th <- theme_pew_serif(base_size = BS) +
  theme(axis.text = element_text(colour = pew_ink), axis.title = element_text(colour = pew_ink),
        strip.text = element_text(face = "bold", colour = pew_ink, size = BS),
        plot.title = element_text(size = BS * 1.15), legend.text = element_text(colour = pew_ink))
short <- function(x) x %>% sub(" \\(.*\\)$", "", .) %>% sub("Hydroxyglutaric acid", "2-Hydroxyglutarate", .) %>%
  sub("Nicotinamide adenine dinucleotide", "NAD", .) %>% sub("3-Hydroxy-3-methylglutaric acid", "HMG", .) %>%
  sub("5-Deoxy-5-\\(methylthio\\)adenosine", "MTA", .)
sym <- function(adj, p) case_when(adj < .05 ~ "*", p < .05 ~ "#", TRUE ~ "")

# ---------- A: heatmap, everything ----------
cellsA <- c(Geno_F3 = "F 3 mo", Geno_F12 = "F 12 mo", Geno_M3 = "M 3 mo", Geno_M12 = "M 12 mo", GxA_F = "F", GxA_M = "M")
get <- function(norm, cons) readRDS(sprintf("results/%s/fit.rds", norm))$wide %>%
  select(Metabolite, matches("_(logFC|P.Value|adj.P.Val)$")) %>%
  pivot_longer(-Metabolite, names_to = c("contrast", ".value"), names_pattern = "(.*)_(logFC|P.Value|adj.P.Val)") %>%
  filter(contrast %in% cons)
dA <- bind_rows(get(CMP, names(cellsA)[1:4]) %>% mutate(block = paste0(nlab[CMP], "\n3xTg vs WT")),
                get(NORM, names(cellsA)[1:4]) %>% mutate(block = paste0(nlab[NORM], "\n3xTg vs WT")),
                get(NORM, names(cellsA)[5:6]) %>% mutate(block = paste0(nlab[NORM], "\nGenotype x age"))) %>%
  mutate(block = factor(block, unique(block)), lab = factor(cellsA[contrast], unique(cellsA)), s = sym(adj.P.Val, P.Value))
m <- dA %>% mutate(k = paste(block, contrast)) %>% select(Metabolite, k, logFC) %>%
  pivot_wider(names_from = k, values_from = logFC) %>% as.data.frame()
rownames(m) <- m$Metabolite; m <- as.matrix(m[, -1])
ord <- rownames(m)[hclust(dist(m), "ward.D2")$order]
dA <- dA %>% mutate(Metabolite = factor(short(Metabolite), unique(short(ord))))
limA <- 1.2
pA <- ggplot(dA, aes(lab, Metabolite, fill = pmax(pmin(logFC, limA), -limA))) +
  geom_tile(colour = "white", linewidth = .5) +
  geom_text(aes(label = s), size = gs * .85, colour = pew_ink, vjust = .75) +
  facet_grid(~block, scales = "free_x", space = "free_x") +
  scale_fill_gradient2("log2 FC\n(3xTg / WT)", low = dn, mid = mid, high = up, limits = c(-limA, limA)) +
  scale_x_discrete(position = "top") +
  labs(title = "A   3xTg vs WT, all 92 metabolites", x = NULL, y = NULL) + th +
  theme(axis.text.y = element_text(size = BS * .7, colour = pew_ink), strip.placement = "outside",
        legend.position = "bottom", legend.key.width = unit(2.2, "lines"), legend.title = element_text(vjust = .8))

# ---------- B: up / down counts, females ----------
cellsF <- c(Geno_F3 = "3 mo", Geno_F12 = "12 mo")
lf <- f$res %>% filter(contrast %in% names(cellsF)) %>% mutate(age = factor(cellsF[contrast], rev(cellsF)))
cnt <- lf %>% group_by(age) %>%
  summarise(up_p = sum(P.Value < .05 & logFC > 0), dn_p = -sum(P.Value < .05 & logFC < 0),
            up_q = sum(adj.P.Val < .05 & logFC > 0), dn_q = -sum(adj.P.Val < .05 & logFC < 0)) %>%
  pivot_longer(-age) %>% mutate(dir = ifelse(grepl("^up", name), "Up in 3xTg", "Down in 3xTg"), lvl = sub("^.._", "", name))
pB <- ggplot(filter(cnt, lvl == "p"), aes(value, age, fill = dir)) +
  geom_col(width = .62, alpha = .35) +
  geom_col(data = filter(cnt, lvl == "q"), width = .62) +
  geom_vline(xintercept = 0, colour = pew_ink) +
  geom_text(data = filter(cnt, lvl == "p", value != 0), aes(label = abs(value), hjust = ifelse(value > 0, -.3, 1.3)),
            size = gs, colour = pew_ink) +
  scale_fill_manual(NULL, values = c(`Up in 3xTg` = up, `Down in 3xTg` = dn), breaks = c("Down in 3xTg", "Up in 3xTg")) +
  scale_x_continuous(labels = abs, expand = expansion(.15)) +
  labs(title = "B   Metabolites changed in 3xTg females", x = "metabolites  (← down  |  up →)", y = NULL) + th +
  theme(legend.position = "top", legend.justification = "left")

# ---------- C: ORA dot plot, female columns ----------
colsC <- c(Geno_F3 = "3 mo", Geno_F12 = "12 mo", GxA_F = "Genotype\nx age")
dbl <- c(KEGG = "KEGG", Reactome_lowest = "Reactome", SMPDB_metabolic = "SMPDB", RefMet_main_class = "RefMet", RefMet_sub_class = "RefMet")
allc <- c("Geno_main", "Geno_F3", "Geno_F12", "Geno_M3", "Geno_M12", "GxA_F", "GxA_M")
er <- read.csv(sprintf("results/enrichment_multi_db_%s.csv", NORM)) %>% filter(contrast %in% allc) %>%
  group_by(database, id) %>% mutate(p = pmin(ora_up_p, ora_down_p), fdr = ifelse(ora_up_p <= ora_down_p, ora_up_fdr, ora_down_fdr)) %>%
  filter(any(p < .05)) %>% ungroup() %>% filter(contrast %in% names(colsC)) %>%
  mutate( db = factor(dbl[database], unique(dbl)), name = sub(" - Mus musculus.*$", "", name),
         lab = sprintf("%s (%d)", ifelse(nchar(name) > 40, paste0(substr(name, 1, 38), "…"), name), n),
         col = factor(colsC[contrast], colsC), s = sym(fdr, p))
ordC <- er %>% group_by(db, lab) %>% summarise(m = min(p), .groups = "drop") %>% arrange(db, desc(m)) %>% pull(lab) %>% unique()
er <- er %>% mutate(lab = factor(lab, ordC), fc = pmax(pmin(mean_logFC, .8), -.8))
pC <- ggplot(er, aes(col, lab)) +
  geom_point(data = filter(er, p >= .05), aes(size = -log10(p), fill = fc), shape = 21, colour = "white", stroke = .3) +
  geom_point(data = filter(er, p < .05), aes(size = -log10(p), fill = fc), shape = 21, colour = pew_ink, stroke = 1.4) +
  geom_text(aes(label = s), size = gs, colour = pew_ink, nudge_x = .3, vjust = .4) +
  facet_grid(db ~ ., scales = "free_y", space = "free_y", switch = "y") +
  scale_fill_gradient2("Mean log2 FC\nof members", low = dn, mid = mid, high = up, limits = c(-.8, .8)) +
  scale_size_continuous("-log10 p\n(ORA)", range = c(2, 9), limits = c(0, 2.2), breaks = c(.5, 1, 2)) +
  guides(size = guide_legend(override.aes = list(colour = pew_ink, fill = mid)), fill = guide_colourbar(order = 1)) +
  labs(title = "C   Pathway enrichment (ORA), females", x = NULL, y = NULL) + th +
  theme(strip.placement = "outside", strip.text.y.left = element_text(angle = 0, hjust = 1, vjust = 1), panel.spacing = unit(.6, "lines"))

# ---------- D: lysine / 2-HG metabolites, females ----------
metsD <- c(`2-Hydroxyglutarate` = "Hydroxyglutaric acid", Pipecolate = "Pipecolic acid",
           `α-Aminoadipate` = "2-Aminoadipic acid", Glutarylcarnitine = "Glutarylcarnitine")
shortD <- setNames(names(metsD), metsD)
gcol <- c(WT = "#2a78d6", `3xTg` = "#e08214")
dD <- as.data.frame(t(f$Y[metsD, , drop = FALSE])) %>% bind_cols(f$meta) %>% filter(Sex == "F") %>%
  pivot_longer(all_of(unname(metsD)), names_to = "Metabolite", values_to = "y") %>%
  mutate(y = log10(2^y + 1), Genotype = factor(ifelse(Genotype == "Tg", "3xTg", "WT"), c("WT", "3xTg")),
         Age = factor(sub("mo", "", Age), c("3", "12")), lab = factor(shortD[Metabolite], names(metsD)))
sD <- dD %>% group_by(lab, Age, Genotype) %>% summarise(m = mean(y), se = sd(y) / sqrt(n()), .groups = "drop")
stD <- f$res %>% filter(Metabolite %in% metsD, contrast %in% names(cellsF)) %>%
  mutate(Age = factor(ifelse(contrast == "Geno_F3", "3", "12"), c("3", "12")), lab = factor(shortD[Metabolite], names(metsD)),
         s = sym(adj.P.Val, P.Value)) %>%
  left_join(dD %>% group_by(lab) %>% summarise(top = max(y)), by = "lab")
pdg <- position_dodge(.35)
pD <- ggplot(sD, aes(Age, m, colour = Genotype, group = Genotype)) +
  geom_point(data = dD, aes(y = y), position = position_jitterdodge(.1, dodge.width = .35, seed = 1), size = 1.8, alpha = .4) +
  geom_line(position = pdg, linewidth = 1.2) +
  geom_errorbar(aes(ymin = m - se, ymax = m + se), width = 0, position = pdg, linewidth = 1) +
  geom_point(position = pdg, size = 3.2) +
  geom_text(data = stD, aes(Age, top, label = s), inherit.aes = FALSE, vjust = -.3, size = gs * 1.3, colour = pew_ink) +
  facet_wrap(~lab, nrow = 1, scales = "free_y") +
  scale_colour_manual(NULL, values = gcol) +
  scale_y_continuous(expand = expansion(c(.08, .25))) +
  labs(title = "D   Lysine catabolism and 2-HG, females", x = "Age (months)", y = expression(log[10]~"(abundance + 1)")) + th +
  theme(legend.position = "top", legend.justification = "left", legend.key.width = unit(2, "lines"),
        panel.border = element_rect(colour = pew_ink, fill = NA, linewidth = .8),
        panel.spacing.x = unit(1.2, "lines"))

# ---------- E: lysine degradation pathway, F 3 mo and F 12 mo ----------
nodes <- tribble(
  ~id,    ~label,               ~x,   ~y,   ~Metabolite,
  "lys",  "Lysine",             3.0,  7.0,  NA,
  "sac",  "Saccharopine",       1.0,  5.6,  NA,
  "pip",  "Pipecolate",         5.0,  5.6,  "Pipecolic acid",
  "aasa", "AASA / P6C",         3.0,  4.2,  NA,
  "aaa",  "α-Aminoadipate",     3.0,  3.0,  "2-Aminoadipic acid",
  "oxo",  "2-Oxoadipate",       3.0,  1.8,  NA,
  "gcoa", "Glutaryl-CoA",       3.0,  0.6,  NA,
  "gcar", "Glutarylcarnitine",  6.0,  0.6,  "Glutarylcarnitine",
  "acoa", "Acetyl-CoA → TCA",   3.0, -0.6,  NA,
  "akg",  "α-Ketoglutarate",   -0.6,  7.0,  NA,
  "hg",   "2-Hydroxyglutarate",-0.6,  8.3,  "Hydroxyglutaric acid")
edges <- tribble(~from, ~to, ~type,
  "lys", "sac", "main", "lys", "pip", "main", "sac", "aasa", "main", "pip", "aasa", "main",
  "aasa", "aaa", "main", "aaa", "oxo", "main", "oxo", "gcoa", "main", "gcoa", "acoa", "main",
  "gcoa", "gcar", "side", "akg", "hg", "side") %>%
  left_join(select(nodes, from = id, x, y), by = "from") %>% left_join(select(nodes, to = id, xend = x, yend = y), by = "to") %>%
  mutate(dx = xend - x, dy = yend - y, L = sqrt(dx^2 + dy^2), sx = ifelse(abs(dx) > abs(dy), .8, .3),
         x = x + dx / L * pmin(sx, L * .3), y = y + dy / L * pmin(.34, L * .3),
         xend = xend - dx / L * pmin(sx, L * .3), yend = yend - dy / L * pmin(.34, L * .3))
nF <- table(f$meta$cell)
pathway <- function(k) {
  st <- f$res %>% filter(contrast == k) %>% select(Metabolite, logFC, P.Value, adj.P.Val)
  nd <- nodes %>% left_join(st, by = "Metabolite") %>%
    mutate(meas = !is.na(logFC), s = ifelse(meas, sym(adj.P.Val, P.Value), ""),
           txt = ifelse(meas, sprintf("%s\n%+.2f%s", label, logFC, s), label),
           fc = logFC, tcol = ifelse(meas & abs(fc) > .8, "white", pew_ink))
  a <- ifelse(k == "Geno_F3", "3", "12")
  ggplot() +
    geom_segment(data = filter(edges, type == "main"), aes(x, y, xend = xend, yend = yend), colour = pew_muted,
                 linewidth = .7, arrow = arrow(length = unit(.18, "cm"), type = "closed")) +
    geom_segment(data = filter(edges, type == "side"), aes(x, y, xend = xend, yend = yend), colour = pew_muted,
                 linewidth = .6, linetype = 2) +
    annotate("curve", x = .5, y = 7.0, xend = 1.9, yend = 6.4, curvature = -.2, colour = pew_muted,
             linewidth = .6, arrow = arrow(length = unit(.15, "cm"), type = "closed")) +
    annotate("text", x = 1.5, y = 4.75, label = "classical (mitochondrial)", size = gs * .7, colour = pew_ink,
             lineheight = .9, hjust = 1, family = "Libre Franklin", fontface = "italic") +
    annotate("text", x = 4.3, y = 4.75, label = "brain-enriched (peroxisomal)", size = gs * .7, colour = pew_ink,
             lineheight = .9, hjust = 0, family = "Libre Franklin", fontface = "italic") +
    geom_label(data = filter(nd, !meas), aes(x, y, label = txt), family = "Libre Franklin", size = gs * .8, colour = pew_ink,
               fill = "white", label.r = unit(.25, "lines"), label.padding = unit(.35, "lines"), linewidth = .25) +
    geom_label(data = filter(nd, meas), aes(x, y, label = txt, fill = fc, colour = I(tcol)), family = "Libre Franklin",
               fontface = "bold", size = gs * .85, lineheight = .95, label.r = unit(.25, "lines"),
               label.padding = unit(.4, "lines"), linewidth = 0) +
    scale_fill_gradientn(colours = c(dn, mid, up), limits = c(-limA, limA), oob = scales::squish, guide = "none") +  # same scale as A
    scale_x_continuous(limits = c(-1.7, 7.1)) + scale_y_continuous(limits = c(-1, 8.8)) +
    labs(title = sprintf("Female, %s mo  (3xTg n = %d, WT n = %d)", a, nF[paste0("Tg_F_", a, "mo")], nF[paste0("WT_F_", a, "mo")])) +
    coord_cartesian(clip = "off") + theme_void(base_family = "Libre Franklin", base_size = BS) +
    theme(plot.title = element_text(face = "bold", colour = pew_ink, size = BS, hjust = 0),
          legend.title = element_text(colour = pew_ink), legend.text = element_text(colour = pew_ink))
}
pE <- pathway("Geno_F3") | pathway("Geno_F12")

# ---------- assemble ----------
right <- wrap_elements(full = pB) / wrap_elements(full = pC) / wrap_elements(full = pD) /
  wrap_elements(full = pE + plot_annotation(title = "E   Lysine degradation pathway, females",
                                            theme = theme(plot.title = element_text(family = "Source Serif 4", face = "bold",
                                                                                    size = BS * 1.15, colour = pew_ink)))) +
  plot_layout(heights = c(2.4, 5.2, 4.6, 7.4))
g <- (wrap_elements(full = pA) | right) + plot_layout(widths = c(.9, 1.1)) +
  plot_annotation(theme = theme(plot.background = element_rect(fill = "white", colour = NA)))
fn <- "figures/figure6"
ggsave(paste0(fn, ".png"), g, width = 24, height = 26, dpi = 200, bg = "white", device = ragg::agg_png)
ggsave(paste0(fn, ".pdf"), g, width = 24, height = 26, device = cairo_pdf, bg = "white")
cat("done\n")
