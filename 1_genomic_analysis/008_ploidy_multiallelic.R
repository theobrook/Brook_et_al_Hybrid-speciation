################################################################################
# 008_ploidy_multiallelic.R
# Brook et al. (2026) Not a dead end: hybrid speciation and syngameon-driven
# dynamics in an island plant radiation
#
# Ploidy check from multiallelic sites (Methods S6; Results S3; Table S16).
# Multiallelic sites (>= 3 alleles) are expressed as a percentage of all
# callable (non-missing) sites per individual, counted in the merged BCF
# before any filtering (008_ploidy_multiallelic.sh).
#
# Species with diploid chromosome counts (2n = 44; Papadopulos et al.):
#   C. huttoniana, C. lanceolaris, C. ptotopetra, C. savolainenii (n = 33)
# Species of unknown ploidy: C. putida, C. inopinata (n = 13)
#
#   1. Summary statistics per group
#   2. Welch's t-test (unequal variances) and Cohen's d
################################################################################

library(tidyverse)


# ------------------------------------------------------------------------------
# Settings
# ------------------------------------------------------------------------------

multiallelic_file <- "data/multiallelic_sites.csv"   # Table S16
out_dir           <- "results"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

diploid_species <- c("huttoniana", "lanceolaris", "ptotopetra", "savolainenii")


# ------------------------------------------------------------------------------
# 1. Data and summary statistics
# ------------------------------------------------------------------------------

sites <- read_csv(multiallelic_file, show_col_types = FALSE) %>%
  mutate(percent = 100 * multiallelic_sites / total_called_sites,
         ploidy  = factor(if_else(species %in% diploid_species, "diploid", "unknown"),
                          levels = c("diploid", "unknown")))

sites %>%
  group_by(ploidy) %>%
  summarise(n = n(), mean = mean(percent), sd = sd(percent),
            min = min(percent), max = max(percent)) %>%
  write_csv(file.path(out_dir, "multiallelic_summary.csv"))


# ------------------------------------------------------------------------------
# 2. Welch's t-test and effect size
# ------------------------------------------------------------------------------

welch <- t.test(percent ~ ploidy, data = sites, var.equal = FALSE)

x <- sites$percent[sites$ploidy == "diploid"]
y <- sites$percent[sites$ploidy == "unknown"]
pooled_sd <- sqrt(((length(x) - 1) * var(x) + (length(y) - 1) * var(y)) /
                    (length(x) + length(y) - 2))

tibble(t = unname(welch$statistic), df = unname(welch$parameter), p = welch$p.value,
       cohens_d = (mean(x) - mean(y)) / pooled_sd) %>%
  write_csv(file.path(out_dir, "multiallelic_welch_test.csv"))
