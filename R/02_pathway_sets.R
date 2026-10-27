# Build metabolite sets (KEGG, Reactome, SMPDB) for the 92 measured metabolites.
# OPTIONAL: the output annotation/pathway_sets_all_dbs.rds is shipped with the repo. Re-running needs the
# third-party downloads listed in README.md placed in annotation/db/ (not version-controlled, ~400 MB).
suppressPackageStartupMessages({library(dplyr); library(tidyr)})
DB <- "annotation/db"
mets <- rownames(readRDS("results/pqn/fit.rds")$Y)
idmap <- read.csv("annotation/kegg_id_map.csv", check.names = FALSE) %>% separate_rows(KEGG, sep = ";") %>% filter(KEGG != "")

# --- ChEBI IDs: KEGG conv, expanded over conjugate acid/base + tautomer links (2 hops) ---
k2c <- read.delim(file.path(DB, "kegg2chebi.tsv"), header = FALSE, col.names = c("KEGG", "chebi")) %>%
  mutate(KEGG = sub("cpd:", "", KEGG), chebi = as.integer(sub("chebi:", "", chebi)))
rel <- read.delim(gzfile(file.path(DB, "relation.tsv.gz"))) %>% filter(relation_type_id %in% c(6, 7, 11)) %>% select(a = init_id, b = final_id)
rel <- bind_rows(rel, rename(rel, a = b, b = a))
chebi <- idmap %>% inner_join(k2c, by = "KEGG", relationship = "many-to-many") %>% distinct(Metabolite, chebi)
for (i in 1:2) chebi <- bind_rows(chebi, chebi %>% inner_join(rel, by = c("chebi" = "a"), relationship = "many-to-many") %>% transmute(Metabolite, chebi = b)) %>% distinct()
write.csv(chebi, "annotation/chebi_id_map_expanded.csv", row.names = FALSE)

sets <- list()
# KEGG (mouse, no overview maps)
kg <- readRDS("annotation/kegg_mmu_pathway_compounds.rds")
sets$KEGG <- kg$link %>% filter(map %in% kg$mmu$map, !grepl("^map01[1-2]", map)) %>% inner_join(idmap, by = c("cpd" = "KEGG"), relationship = "many-to-many") %>%
  left_join(kg$mmu, by = "map") %>% transmute(id = map, name, Metabolite)
# Reactome, mouse: lowest-level and all-levels (includes parent pathways)
rd <- function(f) read.delim(file.path(DB, f), header = FALSE, quote = "", col.names = c("chebi", "id", "url", "name", "evidence", "species")) %>%
  filter(species == "Mus musculus") %>% inner_join(chebi, by = "chebi", relationship = "many-to-many") %>% distinct(id, name, Metabolite)
sets$Reactome_lowest <- rd("ChEBI2Reactome.txt"); sets$Reactome_all_levels <- rd("ChEBI2Reactome_All_Levels.txt")
# SMPDB metabolic pathways (human; drop lipid-species-specific pathways)
sm <- bind_rows(lapply(list.files(file.path(DB, "smpdb"), full.names = TRUE), function(f)
  read.csv(f, colClasses = "character", check.names = FALSE)[, c("SMPDB ID", "Pathway Name", "Pathway Subject", "KEGG ID", "ChEBI ID")]))
names(sm) <- c("id", "name", "subject", "KEGG", "chebi")
sm <- sm %>% filter(subject == "Metabolic", !grepl("[0-9]+:[0-9]+", name))
sets$SMPDB_metabolic <- bind_rows(sm %>% select(-chebi) %>% inner_join(idmap, by = "KEGG", relationship = "many-to-many"),
  sm %>% mutate(chebi = suppressWarnings(as.integer(chebi))) %>% inner_join(chebi, by = "chebi", relationship = "many-to-many")) %>% distinct(id, name, Metabolite)
saveRDS(sets, "annotation/pathway_sets_all_dbs.rds")

cov <- bind_rows(lapply(names(sets), function(d) { s <- sets[[d]] %>% group_by(id) %>% summarise(n = n_distinct(Metabolite))
  tibble(database = d, metabolites_covered = n_distinct(sets[[d]]$Metabolite), pathways_total = nrow(s),
         pathways_ge3 = sum(s$n >= 3), pathways_ge5 = sum(s$n >= 5), median_measured_per_pathway_ge3 = median(s$n[s$n >= 3])) }))
print(as.data.frame(cov))
write.csv(cov, "results/database_coverage.csv", row.names = FALSE)
cat("\nNot covered by any pathway DB:", paste(setdiff(mets, unlist(lapply(sets, `[[`, "Metabolite"))), collapse = "; "), "\n")
