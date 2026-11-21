################################################################################
# 004_fst_niche_correlations.R
# Brook et al. (2026) Not a dead end: hybrid speciation and syngameon-driven
# dynamics in an island plant radiation
#
# Tests the syngameon expectation that ecological differentiation and RI are
# not positively correlated (Methods S3, "Testing the relationship between
# ecological differentiation and RI").
#
# Pearson correlations between phylogenetically corrected pairwise FST
# (FST / patristic distance; 1_genomic_analysis/007_fst_hybrid_test.R) and:
#   - PCA centroid distance (002_environmental_pca.R)
#   - Schoener's D and Warren's I niche overlap (003_niche_overlap.R)
################################################################################

library(tidyverse)


# ------------------------------------------------------------------------------
# Settings
# ------------------------------------------------------------------------------

fst_file      <- "../1_genomic_analysis/results/fst_corrected_pairwise.csv"
centroid_file <- "results/pca_centroid_distance_matrix.csv"
schoener_file <- "results/schoeners_D_matrix.csv"
warren_file   <- "results/warrens_I_matrix.csv"
out_dir       <- "results"


# ------------------------------------------------------------------------------
# Data: one row per species pair
# ------------------------------------------------------------------------------

read_matrix <- function(path) as.matrix(read.csv(path, row.names = 1, check.names = FALSE))

centroid <- read_matrix(centroid_file)
schoener <- read_matrix(schoener_file)
warren   <- read_matrix(warren_file)

pairs <- read_csv(fst_file, show_col_types = FALSE) %>%
  rowwise() %>%
  mutate(centroid_distance = centroid[species1, species2],
         schoeners_D       = schoener[species1, species2],
         warrens_I         = warren[species1, species2]) %>%
  ungroup()

write_csv(pairs, file.path(out_dir, "fst_niche_pairwise.csv"))


# ------------------------------------------------------------------------------
# Pearson correlations
# ------------------------------------------------------------------------------

niche_metrics <- c("centroid_distance", "schoeners_D", "warrens_I")

correlations <- map_dfr(niche_metrics, function(m) {
  ct <- cor.test(pairs$fst_corrected, pairs[[m]], method = "pearson")
  tibble(niche_metric = m,
         r            = unname(ct$estimate),
         p_value      = ct$p.value,
         n_pairs      = nrow(pairs))
})

write_csv(correlations, file.path(out_dir, "fst_niche_correlations.csv"))
