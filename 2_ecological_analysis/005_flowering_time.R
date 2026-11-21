################################################################################
# 005_flowering_time.R
# Brook et al. (2026) Not a dead end: hybrid speciation and syngameon-driven
# dynamics in an island plant radiation
#
# Differences in flowering time among species (Methods S5; Results S2;
# Table S15; Figs S8 and S9).
#
#   1. Flowering records from herbarium specimens (NSW and K), field
#      collections and field observations, with a known month and day; each
#      record converted to the number of days after 1 April (flowering season
#      runs from April, so January-March records belong to the following year)
#   2. Kruskal-Wallis rank sum test, followed by pairwise Wilcoxon rank-sum
#      tests with Bonferroni correction across the 15 species pairs:
#        a. all flowering records (Table S15)
#        b. presence-only data, one record per species per day
#   3. Flowering-time densities and medians (Fig. S8)
#   4. Daily presence/absence of flowering from 1 June (Fig. S9)
################################################################################

library(tidyverse)


# ------------------------------------------------------------------------------
# Settings
# ------------------------------------------------------------------------------

flowering_file <- "data/flowering_time.csv"
out_dir        <- "results/005_flowering_time"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

species_cols <- c("C. huttoniana" = "#66C2A5", "C. inopinata" = "#FFD92F",
                  "C. lanceolaris" = "#FC8D62", "C. ptotopetra" = "#A6D854",
                  "C. putida" = "#8DA0CB", "C. savolainenii" = "#E78AC3")

# Season start for each analysis: dates before this month belong to the
# following year
season_date <- function(month_num, day, start_month) {
  as.Date(paste(ifelse(month_num >= start_month, 2026, 2027), month_num, day, sep = "-"))
}


# ------------------------------------------------------------------------------
# 1. Flowering records with a known month and day
# ------------------------------------------------------------------------------

flowering <- read_csv(flowering_file, col_types = cols(.default = "c")) %>%
  filter(Flowers == "Yes", !Month %in% c(NA, "Uncertain"), !Day %in% c(NA, "Uncertain")) %>%
  mutate(month_num = match(Month, month.name),
         Day       = as.integer(Day),
         date      = season_date(month_num, Day, start_month = 4),
         days_after_1_April = as.numeric(date - as.Date("2026-04-01")))

flowering %>%
  group_by(Species) %>%
  summarise(n = n(), median_days_after_1_April = median(days_after_1_April)) %>%
  write_csv(file.path(out_dir, "flowering_summary.csv"))


# ------------------------------------------------------------------------------
# 2. Kruskal-Wallis and pairwise Wilcoxon rank-sum tests
# ------------------------------------------------------------------------------

flowering_tests <- function(d) {
  kw <- kruskal.test(days_after_1_April ~ Species, data = d)

  pairwise <- combn(sort(unique(d$Species)), 2, simplify = FALSE) %>%
    map_dfr(function(p) {
      x <- d$days_after_1_April[d$Species == p[1]]
      y <- d$days_after_1_April[d$Species == p[2]]
      w <- wilcox.test(x, y, exact = FALSE)
      tibble(species_1 = p[1], species_2 = p[2], n1 = length(x), n2 = length(y),
             W = unname(w$statistic), p_raw = w$p.value)
    }) %>%
    mutate(p_adj = p.adjust(p_raw, method = "bonferroni"))

  list(kruskal = tibble(chi_squared = unname(kw$statistic), df = unname(kw$parameter),
                        p = kw$p.value),
       pairwise = pairwise)
}

records  <- flowering_tests(flowering)
presence <- flowering_tests(distinct(flowering, Species, date, .keep_all = TRUE))

write_csv(records$kruskal,   file.path(out_dir, "kruskal_records.csv"))
write_csv(records$pairwise,  file.path(out_dir, "pairwise_records.csv"))     # Table S15
write_csv(presence$kruskal,  file.path(out_dir, "kruskal_presence_only.csv"))
write_csv(presence$pairwise, file.path(out_dir, "pairwise_presence_only.csv"))


# ------------------------------------------------------------------------------
# 3. Flowering-time densities and medians, days after 1 April (Fig. S8)
# ------------------------------------------------------------------------------
# Gaussian kernel, bandwidth 15 days

density_df <- flowering %>%
  group_by(Species) %>%
  reframe({
    d <- density(days_after_1_April, bw = 15, n = 512, from = 0, to = 280)
    tibble(x = d$x, y = d$y)
  })

medians <- flowering %>%
  group_by(Species) %>%
  summarise(median_day = median(days_after_1_April)) %>%
  left_join(density_df %>% group_by(Species) %>% summarise(fx = list(approxfun(x, y))),
            by = "Species") %>%
  mutate(height = map2_dbl(fx, median_day, ~ .x(.y)))

fig_s8 <- ggplot() +
  geom_line(data = density_df, aes(x, y, colour = Species)) +
  geom_segment(data = medians, aes(x = median_day, xend = median_day, y = 0, yend = height,
                                   colour = Species), linetype = "dashed", linewidth = 0.7) +
  scale_colour_manual(values = species_cols) +
  scale_x_continuous(breaks = seq(0, 280, 30)) +
  labs(x = "Number of days after 1 April", y = "Density of flowering events", colour = "Species") +
  theme_minimal(base_size = 12) +
  theme(legend.text = element_text(face = "italic"))

ggsave(file.path(out_dir, "fig_s8_flowering_density.png"), fig_s8,
       width = 15, height = 6, dpi = 300, bg = "white")


# ------------------------------------------------------------------------------
# 4. Daily presence/absence of flowering from 1 June (Fig. S9)
# ------------------------------------------------------------------------------

daily <- flowering %>%
  mutate(day_index = as.numeric(season_date(month_num, Day, start_month = 6) -
                                  as.Date("2026-06-01"))) %>%
  distinct(Species, day_index)

grid <- expand_grid(Species = names(species_cols),
                    day_index = 0:max(daily$day_index)) %>%
  filter(Species %in% daily$Species) %>%
  left_join(mutate(daily, flowering = TRUE), by = c("Species", "day_index")) %>%
  mutate(flowering = replace_na(flowering, FALSE),
         Species   = factor(Species, levels = rev(names(species_cols))))

fig_s9 <- ggplot(grid, aes(day_index, Species)) +
  geom_tile(aes(fill = Species, alpha = flowering), colour = "grey85", linewidth = 0.2) +
  scale_fill_manual(values = species_cols, guide = "none") +
  scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = 0), guide = "none") +
  scale_x_continuous(expand = c(0, 0)) +
  labs(x = "Days after 1 June", y = NULL) +
  theme_minimal(base_size = 12) +
  theme(axis.text.y = element_text(face = "italic"), panel.grid = element_blank())

ggsave(file.path(out_dir, "fig_s9_flowering_presence.png"), fig_s9,
       width = 12, height = 2.5, dpi = 300, bg = "white")
