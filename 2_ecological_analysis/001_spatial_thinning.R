################################################################################
# 001_spatial_thinning.R
# Brook et al. (2026) Not a dead end: hybrid speciation and syngameon-driven
# dynamics in an island plant radiation
#
# Spatially thinned occurrence dataset for the ecological differentiation
# analysis (Methods S3, "Environmental layers and species occurrence data";
# Table S18). One record per species is retained at random within each
# 10 x 10 m grid cell, matching the resolution of the environmental layers.
# The species distribution models used the full (unthinned) dataset.
#
# Coordinates are projected (metres), so thinning is grid-based.
################################################################################

library(tidyverse)


# ------------------------------------------------------------------------------
# Settings
# ------------------------------------------------------------------------------

occ_file <- "data/species_occurrence.csv"
out_dir  <- "species_occurrence/spatial_thinning_10m"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

x_col     <- "Ddlongitude"   # easting (m)
y_col     <- "Ddlatitude"    # northing (m)
cell_size <- 10              # m

set.seed(42)


# ------------------------------------------------------------------------------
# Thinning: one random record per species per grid cell
# ------------------------------------------------------------------------------

occ <- read_csv(occ_file, show_col_types = FALSE) %>%
  drop_na(all_of(c(x_col, y_col)))

# Species are thinned one at a time, in the order they appear in the file, so
# that the random draws reproduce the dataset used in the paper
thinned <- map_dfr(unique(occ$Species), function(sp) {
  occ %>%
    filter(Species == sp) %>%
    mutate(cell = paste(floor(.data[[x_col]] / cell_size),
                        floor(.data[[y_col]] / cell_size), sep = "_")) %>%
    group_by(cell) %>%
    slice_sample(n = 1) %>%
    ungroup() %>%
    select(-cell)
})

write_csv(thinned, file.path(out_dir, "occurrences_thinned.csv"))


# ------------------------------------------------------------------------------
# Records per species before and after thinning (Table S18)
# ------------------------------------------------------------------------------

count(occ, Species, name = "full_dataset") %>%
  left_join(count(thinned, Species, name = "thinned_dataset"), by = "Species") %>%
  write_csv(file.path(out_dir, "thinning_summary.csv"))
