################################################################################
# 007_pixy_summary.R
# Brook et al. (2026) Not a dead end: hybrid speciation and syngameon-driven
# dynamics in an island plant radiation
#
# Genome-wide nucleotide diversity (pi) per species and absolute divergence
# (DXY) per species pair from per-site pixy output (Methods S2,
# "Genome-wide population genetic statistics"; Tables S2 and S3).
# Run after steps 4 and 5 of 007_genetic_diversity.sh.
#
# Genome-wide values are calculated as sum(count_diffs) / sum(count_comparisons)
# across all sites, rather than by averaging per-site values.
################################################################################

library(tidyverse)


# ------------------------------------------------------------------------------
# Settings
# ------------------------------------------------------------------------------

pi_dir  <- "pi"      # pixy pi output
dxy_dir <- "dxy"     # pixy DXY output
out_dir <- "results"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)


# ------------------------------------------------------------------------------
# Function: combine pixy files and sum differences / comparisons by group
# ------------------------------------------------------------------------------

summarise_pixy <- function(dir, pattern, group_cols, stat_name) {
  list.files(dir, pattern = pattern, recursive = TRUE, full.names = TRUE) %>%
    map_dfr(read_tsv, show_col_types = FALSE) %>%
    filter(!is.na(count_diffs), !is.na(count_comparisons), count_comparisons > 0) %>%
    group_by(across(all_of(group_cols))) %>%
    summarise(n_sites         = n(),
              sum_diffs       = sum(count_diffs),
              sum_comparisons = sum(count_comparisons),
              "{stat_name}"   := sum_diffs / sum_comparisons,
              .groups = "drop")
}


# ------------------------------------------------------------------------------
# Nucleotide diversity per species
# ------------------------------------------------------------------------------

pi_per_species <- summarise_pixy(pi_dir, "_pi\\.txt$", "pop", "pi")
write_tsv(pi_per_species, file.path(out_dir, "pixy_pi_per_species.tsv"))


# ------------------------------------------------------------------------------
# Absolute divergence per species pair
# ------------------------------------------------------------------------------

dxy_per_pair <- summarise_pixy(dxy_dir, "_dxy\\.txt$", c("pop1", "pop2"), "dxy")
write_tsv(dxy_per_pair, file.path(out_dir, "pixy_dxy_per_pair.tsv"))
