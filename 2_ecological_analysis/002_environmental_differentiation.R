################################################################################
# 002_environmental_differentiation.R
# Brook et al. (2026) Not a dead end: hybrid speciation and syngameon-driven
# dynamics in an island plant radiation
#
# Tests of environmental differentiation among species (Methods S3,
# "Ecological differentiation analysis"; Tables S7 and S8).
# Run after 002_environmental_pca.R.
#
#   1. Global PERMANOVA ~ species on the standardised environmental variables
#      (adonis2; Euclidean distance; 9,999 permutations)
#   2. Pairwise PERMANOVAs for the 15 species pairs, Bonferroni-corrected, with
#      the Euclidean distance between species centroids across all nine PCs
#      (Table S7)
#   3. Pairwise Wilcoxon rank-sum tests on the scores of each of the nine PCs,
#      Bonferroni-corrected across 135 tests (15 pairs x 9 PCs; Table S8)
#
# Software: R v.4.x, vegan v.2.7.1
################################################################################

library(tidyverse)
library(vegan)


# ------------------------------------------------------------------------------
# Settings
# ------------------------------------------------------------------------------

scores_file <- "results/data_with_pc_scores.csv"   # from 002_environmental_pca.R
out_dir     <- "results"

env_vars <- c("aspect", "coastdist", "creekdist", "dem", "geology",
              "hillshade", "slope", "vegetation", "vrm")
pc_cols  <- paste0("PC", 1:9)
n_perm   <- 9999


# ------------------------------------------------------------------------------
# Data: standardised environmental variables and PC scores
# ------------------------------------------------------------------------------

env <- read_csv(scores_file, show_col_types = FALSE) %>%
  mutate(Species = factor(Species))

X     <- scale(env[, env_vars])
pairs <- combn(levels(env$Species), 2, simplify = FALSE)


# ------------------------------------------------------------------------------
# 1. Global PERMANOVA
# ------------------------------------------------------------------------------

set.seed(42)
global <- adonis2(X ~ Species, data = env, method = "euclidean", permutations = n_perm)

tibble(F = global$F[1], R2 = global$R2[1], p = global$`Pr(>F)`[1]) %>%
  write_csv(file.path(out_dir, "permanova_global.csv"))


# ------------------------------------------------------------------------------
# 2. Pairwise PERMANOVAs and centroid distances (Table S7)
# ------------------------------------------------------------------------------

centroids <- env %>%
  group_by(Species) %>%
  summarise(across(all_of(pc_cols), mean)) %>%
  column_to_rownames("Species")

set.seed(42)
pairwise_permanova <- map_dfr(pairs, function(p) {
  idx <- env$Species %in% p
  res <- adonis2(X[idx, ] ~ droplevels(env$Species[idx]),
                 method = "euclidean", permutations = n_perm)
  tibble(species_1         = p[1],
         species_2         = p[2],
         centroid_distance = sqrt(sum((centroids[p[1], ] - centroids[p[2], ])^2)),
         F                 = res$F[1],
         R2                = res$R2[1],
         p_perm            = res$`Pr(>F)`[1])
}) %>%
  mutate(p_adj = p.adjust(p_perm, method = "bonferroni"))

write_csv(pairwise_permanova, file.path(out_dir, "permanova_pairwise.csv"))


# ------------------------------------------------------------------------------
# 3. Pairwise Wilcoxon tests on each PC (Table S8)
# ------------------------------------------------------------------------------

pc_wilcoxon <- expand_grid(pc = pc_cols, pair = pairs) %>%
  pmap_dfr(function(pc, pair) {
    x <- env[[pc]][env$Species == pair[1]]
    y <- env[[pc]][env$Species == pair[2]]
    w <- wilcox.test(x, y)
    tibble(pc = pc, species_1 = pair[1], species_2 = pair[2],
           W = unname(w$statistic), p_raw = w$p.value,
           mean_diff = mean(x) - mean(y))
  }) %>%
  mutate(p_adj = p.adjust(p_raw, method = "bonferroni"))

write_csv(pc_wilcoxon, file.path(out_dir, "pc_wilcoxon_pairwise.csv"))
