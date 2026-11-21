################################################################################
# 002_environmental_pca.R
# Brook et al. (2026) Not a dead end: hybrid speciation and syngameon-driven
# dynamics in an island plant radiation
#
# PCA of environmental variables at occurrence points and pairwise distances
# between species centroids in PC space (Methods S3, "Ecological
# differentiation analysis").
#
#   1. Standardise the nine environmental variables (zero mean, unit variance)
#   2. PCA (prcomp)
#   3. Species centroids in PC space and pairwise Euclidean distances between
#      centroids across all nine PCs (Table S7)
#   4. PCA figures with 95% confidence ellipses and loading arrows
#      (Fig. 3: PC1 v PC2, PC1 v PC3; Fig. S4: PC1 v PC2, PC3, PC4 and scree plot)
#
# Input: environmental values at the spatially thinned occurrence points, from
# scripts/extract_values.py (run after 001_spatial_thinning.R).
# The centroid distance matrix is used in 004_fst_niche_correlations.R.
################################################################################

library(tidyverse)
library(ggrepel)
library(patchwork)


# ------------------------------------------------------------------------------
# Settings
# ------------------------------------------------------------------------------

env_data_file <- "species_occurrence/spatial_thinning_10m/extracted_data.csv"   # spatially thinned occurrences (Table S18)
out_dir       <- "results"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

focal_species <- c("huttoniana", "inopinata", "lanceolaris",
                   "ptotopetra", "putida", "savolainenii")

env_vars <- c("aspect", "coastdist", "creekdist", "dem", "geology",
              "hillshade", "slope", "vegetation", "vrm")

env_labels <- c(aspect = "Aspect", coastdist = "Distance to coast",
                creekdist = "Distance to creek", dem = "Elevation",
                geology = "Geology", hillshade = "Hillshade",
                slope = "Gradient of slope", vegetation = "Vegetation", vrm = "VRM")

species_colours <- c(huttoniana = "#66C2A5", inopinata = "#FFD92F",
                     lanceolaris = "#FC8D62", ptotopetra = "#A6D854",
                     putida = "#8DA0CB", savolainenii = "#E78AC3")


# ------------------------------------------------------------------------------
# 1. Data
# ------------------------------------------------------------------------------

env <- read_csv(env_data_file, show_col_types = FALSE) %>%
  filter(Species %in% focal_species) %>%
  drop_na(all_of(env_vars))


# ------------------------------------------------------------------------------
# 2. PCA on standardised variables
# ------------------------------------------------------------------------------

pca <- prcomp(scale(env[, env_vars]), center = FALSE, scale. = FALSE)

var_explained <- summary(pca)$importance
prop_var      <- var_explained["Proportion of Variance", ]

env <- bind_cols(env, as_tibble(pca$x))

write_csv(as_tibble(var_explained, rownames = "statistic"),
          file.path(out_dir, "pca_variance_explained.csv"))
write_csv(as_tibble(pca$rotation, rownames = "variable"),
          file.path(out_dir, "pca_loadings.csv"))
write_csv(env, file.path(out_dir, "data_with_pc_scores.csv"))


# ------------------------------------------------------------------------------
# 3. Species centroids and Euclidean centroid distances (all nine PCs)
# ------------------------------------------------------------------------------

pc_cols <- colnames(pca$x)

centroids <- env %>%
  group_by(Species) %>%
  summarise(across(all_of(pc_cols), mean)) %>%
  column_to_rownames("Species")

centroid_distance <- function(sp1, sp2) {
  d <- as.numeric(centroids[sp1, ] - centroids[sp2, ])
  sqrt(sum(d^2))
}

dist_mat <- outer(focal_species, focal_species, Vectorize(centroid_distance))
dimnames(dist_mat) <- list(focal_species, focal_species)

write.csv(dist_mat, file.path(out_dir, "pca_centroid_distance_matrix.csv"))


# ------------------------------------------------------------------------------
# 4. PCA figures (Fig. 3 and Fig. S4)
# ------------------------------------------------------------------------------

pca_plot <- function(pc_x, pc_y, loading_scale = 3) {
  pct <- round(100 * prop_var[c(pc_x, pc_y)], 1)
  arrows <- tibble(label = env_labels[env_vars],
                   x = pca$rotation[, pc_x] * loading_scale,
                   y = pca$rotation[, pc_y] * loading_scale)

  ggplot(env, aes(.data[[pc_x]], .data[[pc_y]], fill = Species)) +
    stat_ellipse(aes(colour = Species), level = 0.95, geom = "polygon",
                 alpha = 0.12, linetype = "dashed", linewidth = 0.8) +
    geom_point(shape = 21, size = 2.5, alpha = 0.6, stroke = 0.3, colour = "white") +
    geom_segment(data = arrows, inherit.aes = FALSE,
                 aes(x = 0, y = 0, xend = x, yend = y),
                 arrow = arrow(length = unit(0.2, "cm")), colour = "#333333", linewidth = 0.6) +
    geom_text_repel(data = arrows, inherit.aes = FALSE,
                    aes(x = x * 1.15, y = y * 1.15, label = label),
                    size = 2.8, fontface = "bold", colour = "#333333",
                    bg.color = "white", bg.r = 0.15, seed = 42) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey70", linewidth = 0.3) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey70", linewidth = 0.3) +
    scale_fill_manual(values = species_colours, labels = function(x) paste0("C. ", x)) +
    scale_colour_manual(values = species_colours, guide = "none") +
    labs(x = paste0(pc_x, " (", pct[1], "%)"), y = paste0(pc_y, " (", pct[2], "%)")) +
    theme_minimal(base_size = 11) +
    theme(panel.grid.minor = element_blank(), legend.text = element_text(face = "italic"))
}

scree <- tibble(PC = seq_along(prop_var),
                variance = 100 * prop_var,
                cumulative = 100 * var_explained["Cumulative Proportion", ]) %>%
  ggplot(aes(PC, variance)) +
  geom_col(fill = "#4a4a4a", alpha = 0.7, width = 0.6) +
  geom_line(aes(y = cumulative), colour = "#E34234", linewidth = 1) +
  geom_point(aes(y = cumulative), colour = "#E34234", size = 2.5) +
  scale_x_continuous(breaks = seq_along(prop_var)) +
  scale_y_continuous(limits = c(0, 105)) +
  labs(x = "PC", y = "Variance (%)") +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank())

fig_3 <- pca_plot("PC1", "PC2") + pca_plot("PC1", "PC3") +
  plot_layout(guides = "collect") + plot_annotation(tag_levels = "A")
ggsave(file.path(out_dir, "fig3_environmental_pca.png"), fig_3,
       width = 14, height = 6, dpi = 300, bg = "white")

fig_s4 <- pca_plot("PC1", "PC2") + pca_plot("PC1", "PC3") + pca_plot("PC1", "PC4") + scree +
  plot_layout(ncol = 2, guides = "collect") + plot_annotation(tag_levels = "A")
ggsave(file.path(out_dir, "figS4_environmental_pca.png"), fig_s4,
       width = 14, height = 12, dpi = 300, bg = "white")
