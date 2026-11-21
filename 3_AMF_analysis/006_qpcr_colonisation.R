################################################################################
# 006_qpcr_colonisation.R
# Brook et al. (2026) Not a dead end: hybrid speciation and syngameon-driven
# dynamics in an island plant radiation
#
# Differences in root AM fungal abundance among host species from qPCR
# (Methods S4, "Quantification of AM fungal abundance with qPCR"; Results S1;
# Table S22).
#
# Response variables:
#   abs_root_amf_ng        absolute root AM fungal DNA (ng, from the standard curve)
#   norm_abs_root_amf      absolute root AM fungal DNA / total DNA extracted
#   rel_colonisation       efficiency-corrected root AM fungal : plant DNA,
#                          E_plant^Ct_plant / E_fungal^Ct_fungal
#   log10_rel_colonisation log10 of the above
#
#   1. Shapiro-Wilk tests within each species and Levene's test across species
#   2. Assumptions were violated for every response, so differences among
#      species were tested with Kruskal-Wallis rank sum tests, followed by
#      Dunn's post hoc tests with Bonferroni correction
################################################################################

library(tidyverse)
library(car)       # leveneTest
library(rstatix)   # shapiro_test, dunn_test


# ------------------------------------------------------------------------------
# Settings
# ------------------------------------------------------------------------------

qpcr_file <- "data/qpcr_results.csv"
out_dir   <- "results/006_qpcr"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

responses <- c(abs_root_amf_ng        = "Absolute root AM fungal DNA (ng)",
               norm_abs_root_amf      = "Normalised absolute root AM fungal DNA (ng/ng)",
               rel_colonisation       = "Relative colonisation",
               log10_rel_colonisation = "log10(relative colonisation)")


# ------------------------------------------------------------------------------
# Data
# ------------------------------------------------------------------------------

qpcr <- read_csv(qpcr_file, show_col_types = FALSE) %>%
  mutate(species = factor(species))

qpcr %>%
  group_by(species) %>%
  summarise(n = n(), across(all_of(names(responses)),
                            list(mean = ~ mean(.x, na.rm = TRUE), sd = ~ sd(.x, na.rm = TRUE)))) %>%
  write_csv(file.path(out_dir, "summary_by_species.csv"))


# ------------------------------------------------------------------------------
# 1-2. Assumption checks, Kruskal-Wallis and Dunn's tests
# ------------------------------------------------------------------------------

results <- imap(responses, function(label, var) {
  d <- qpcr %>% select(species, value = all_of(var)) %>% drop_na()

  list(shapiro = d %>% group_by(species) %>% shapiro_test(value) %>% mutate(response = label),
       levene  = leveneTest(value ~ species, data = d) %>% broom::tidy() %>% mutate(response = label),
       kruskal = kruskal.test(value ~ species, data = d) %>% broom::tidy() %>% mutate(response = label),
       dunn    = d %>% dunn_test(value ~ species, p.adjust.method = "bonferroni") %>%
                   mutate(response = label))
})

for (test in c("shapiro", "levene", "kruskal", "dunn")) {
  map_dfr(results, test) %>%
    write_csv(file.path(out_dir, paste0(test, ".csv")))
}
