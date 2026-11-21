################################################################################
# 007_fst_hybrid_test.R
# Brook et al. (2026) Not a dead end: hybrid speciation and syngameon-driven
# dynamics in an island plant radiation
#
# Tests whether FST between putative hybrid species and their parents differs
# from FST between non-hybrid species pairs, after correcting for phylogenetic
# relatedness (Methods S2, "Genome-wide population genetic statistics").
#
#   1. Pairwise weighted Weir & Cockerham FST from VCFtools (--weir-fst-pop)
#   2. Mean interspecific patristic distance from the ML phylogeny
#      (008_phylogenetics.sh), excluding C67 and C70
#   3. Corrected FST = FST / mean patristic distance
#   4. Wilcoxon rank-sum test: hybrid-parent vs non-hybrid pairs
#
# Hybrid-parent pairs include both candidate second parents of
# C. savolainenii (C. inopinata and C. lanceolaris).
#
# The corrected FST table is also used in
# 2_ecological_analysis/004_fst_niche_correlations.R.
################################################################################

library(ape)
library(tidyverse)


# ------------------------------------------------------------------------------
# Settings
# ------------------------------------------------------------------------------

fst_file    <- "results/pairwise_fst_summary_ingroup.tsv"   # VCFtools weighted FST
tree_file   <- "ML_search.raxml.bestTree"                   # ML tree (008_phylogenetics.sh)
popmap_file <- "samples.txt"                                # sample, species
out_dir     <- "results"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

focal_species   <- c("huttoniana", "inopinata", "lanceolaris",
                     "ptotopetra", "putida", "savolainenii")
exclude_samples <- c("C67", "C70")

hybrid_parent_pairs <- list(
  c("putida",       "ptotopetra"),
  c("putida",       "huttoniana"),
  c("savolainenii", "ptotopetra"),
  c("savolainenii", "inopinata"),
  c("savolainenii", "lanceolaris")
)


# ------------------------------------------------------------------------------
# 1. Pairwise FST (columns: Pop1 Pop2 N_Pop1 N_Pop2 N_SNPs Weighted_Fst Mean_Fst)
# ------------------------------------------------------------------------------

fst <- read_lines(fst_file, skip = 1) %>%
  str_trim() %>%
  discard(~ .x == "") %>%
  str_split("\\s+", simplify = TRUE) %>%
  as_tibble(.name_repair = ~ c("species1", "species2", "n_pop1", "n_pop2",
                               "n_snps", "weighted_fst", "mean_fst")) %>%
  mutate(across(n_pop1:mean_fst, as.numeric))

get_fst <- function(sp1, sp2) {
  fst %>%
    filter((species1 == sp1 & species2 == sp2) | (species1 == sp2 & species2 == sp1)) %>%
    pull(weighted_fst) %>%
    first()
}


# ------------------------------------------------------------------------------
# 2. Mean interspecific patristic distance
# ------------------------------------------------------------------------------

pdist  <- cophenetic.phylo(read.tree(tree_file))
popmap <- read.table(popmap_file, col.names = c("sample_id", "species")) %>%
  filter(!sample_id %in% exclude_samples, species %in% focal_species)

mean_patristic <- function(sp1, sp2) {
  t1 <- intersect(popmap$sample_id[popmap$species == sp1], rownames(pdist))
  t2 <- intersect(popmap$sample_id[popmap$species == sp2], rownames(pdist))
  mean(pdist[t1, t2, drop = FALSE], na.rm = TRUE)
}


# ------------------------------------------------------------------------------
# 3. Corrected FST for all 15 species pairs
# ------------------------------------------------------------------------------

is_hybrid_parent <- function(sp1, sp2) {
  any(map_lgl(hybrid_parent_pairs, ~ setequal(.x, c(sp1, sp2))))
}

fst_corrected <- combn(focal_species, 2, simplify = FALSE) %>%
  map_dfr(~ tibble(species1 = .x[1], species2 = .x[2])) %>%
  rowwise() %>%
  mutate(fst       = get_fst(species1, species2),
         patristic = mean_patristic(species1, species2),
         pair_type = if_else(is_hybrid_parent(species1, species2),
                             "hybrid_parent", "non_hybrid")) %>%
  ungroup() %>%
  mutate(fst_corrected = fst / patristic,
         pair_type     = factor(pair_type, levels = c("non_hybrid", "hybrid_parent")))

write_csv(fst_corrected, file.path(out_dir, "fst_corrected_pairwise.csv"))


# ------------------------------------------------------------------------------
# 4. Wilcoxon rank-sum test
# ------------------------------------------------------------------------------

wilcox_result <- wilcox.test(fst_corrected ~ pair_type, data = fst_corrected,
                             alternative = "two.sided", exact = FALSE)

capture.output(
  fst_corrected %>%
    group_by(pair_type) %>%
    summarise(n = n(), mean = mean(fst_corrected), sd = sd(fst_corrected),
              median = median(fst_corrected)),
  wilcox_result,
  file = file.path(out_dir, "fst_hybrid_test.txt")
)
