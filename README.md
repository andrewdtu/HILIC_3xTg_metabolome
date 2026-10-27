# HILIC_3xTg_metabolome

Simcox Lab. Analysis code for the 3xTg-AD hippocampal metabolomics figure (Figure 6) in
Mill *et al.*, "A framework for robust targeted metabolomic profiling using chromatographic
reproducibility and rational internal standards" (JASMS).

This repository holds only the code used in the paper: the BCA-protein-normalized data, protein + PQN
normalization, limma, pathway enrichment and Figure 6.

## Design
3xTg-AD vs C57BL/6J (B6J) mice, 2 × 2 × 2 design (genotype × sex × age: 3 and 12 months), n = 4–5 per group, 38 mice,
92 targeted HILIC-MS metabolites.

## Running it
From the repository root:

```bash
./run_all.sh
```

| Step | Script | Output |
|---|---|---|
| 1 | `R/01_factorial_limma.R none` | `results/none/` — limma on protein-normalized data (as supplied) |
| 1 | `R/01_factorial_limma.R pqn` | `results/pqn/` — limma on protein + global PQN (primary) |
| 2 | `R/02_pathway_sets.R` *(optional)* | `annotation/pathway_sets_all_dbs.rds` (shipped; rebuild needs `annotation/db/`) |
| 3 | `R/03_enrichment.R pqn` | `results/enrichment_multi_db_pqn.csv` — ORA and camera, KEGG / Reactome / SMPDB / RefMet |
| 4 | `R/04_figure6.R` | `figures/figure6.png`, `figures/figure6.pdf` |

## Methods
- **Input.** Each metabolite is in ng per µg BCA protein, from `data/JASMS metaboanalyst RUN 1 MAIN + run 2.csv`.
- **Normalization.**
  - Zeros (9 of 3,496 values) are replaced by half of that metabolite's minimum positive value.
  - Probabilistic quotient normalization (PQN; Dieterle *et al.*, *Anal Chem* 2006) is then applied. The reference is the per-metabolite median across all samples, and each sample is divided by the median of its ratios to that reference.
  - Values are then log2-transformed.
- **Statistics.**
  - limma uses a cell-means model of the eight genotype × sex × age groups.
  - Contrasts: 3xTg vs B6J within each sex × age, and genotype × age within each sex. Additional contrasts are written to the results tables.
  - Moderated t-tests use empirical Bayes shrinkage. P-values are adjusted by Benjamini–Hochberg within each contrast.
- **Pathway enrichment (Figure 6C).**
  - Directional over-representation analysis: a one-sided hypergeometric test on metabolites with p < 0.05, run separately for increased and decreased metabolites.
  - Databases: KEGG (mouse), Reactome (mouse, lowest level), SMPDB (metabolic) and RefMet classes.
  - Background: the measured metabolites annotated in each resource. Only sets with ≥ 3 measured members are tested.
- **Symbols:** `*` FDR < 0.05, `#` nominal p < 0.05.
- **Figure 6E.** The lysine degradation scheme is adapted from Tan, Chrysopoulou & Rinschen, *Physiol Genomics* 2023;55:579.

## Software
R 4.5.0 with limma 3.64.3, dplyr 1.2.0, tidyr 1.3.2, ggplot2 4.0.2, patchwork 1.3.2, ragg 1.5.2 and scales 1.4.0.
The figure uses the Libre Franklin and Source Serif 4 fonts (Google Fonts), which must be installed system-wide.

## Annotation files (`annotation/`)
- `kegg_id_map.csv` — curated KEGG compound IDs for the 92 metabolites. 2-HG maps to C02630, C03196 and C01087, because the enantiomers are not resolved.
- `refmet_classes.csv` — RefMet main and sub classes.
- `kegg_mmu_pathway_compounds.rds` — cached KEGG REST pathway–compound links (mouse).
- `pathway_sets_all_dbs.rds` — the metabolite sets used by `03_enrichment.R`.

To rebuild `pathway_sets_all_dbs.rds`, put these downloads in `annotation/db/`, which is not version-controlled:
- Reactome `ChEBI2Reactome.txt` and `ChEBI2Reactome_All_Levels.txt` (reactome.org/download-data)
- ChEBI `relation.tsv.gz` (EBI ChEBI flat files)
- `kegg2chebi.tsv` (`https://rest.kegg.jp/conv/chebi/compound`)
- SMPDB per-pathway metabolite CSVs, unzipped into `annotation/db/smpdb/` (smpdb.ca/downloads, "Metabolite names linked to SMPDB pathways")
