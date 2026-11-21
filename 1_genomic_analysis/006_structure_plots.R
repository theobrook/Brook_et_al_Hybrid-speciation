################################################################################
# 006_structure_plots.R
# Brook et al. (2026) Not a dead end: hybrid speciation and syngameon-driven
# dynamics in an island plant radiation
#
# ADMIXTURE and STRUCTURE bar plots for K = 2-6 (Fig. 1D-E; Methods S2,
# "Model-based ancestry estimation"). For each K, the single replicate with
# the lowest cross-validation (CV) error (ADMIXTURE) or the highest lnP(D)
# (STRUCTURE) out of 50 is plotted. Run after 006_population_structure.sh.
################################################################################

library(tidyverse)
library(patchwork)


# ------------------------------------------------------------------------------
# Settings
# ------------------------------------------------------------------------------

admix_dir   <- "admixture/results"     # K{K}/K{K}_rep{R}.Q and K{K}/cv_K{K}_rep{R}.error
struct_dir  <- "structure/results"     # K{K}/K{K}_rep{R}_f
fam_file    <- "snps.fam"              # sample order of the ADMIXTURE / STRUCTURE input
popmap_file <- "samples.txt"           # sample, species
out_dir     <- "results"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

k_values <- 2:6
n_reps   <- 50

species_cols <- c(inopinata = "#FFD92F", lanceolaris = "#FC8D62", huttoniana = "#66C2A5",
                  savolainenii = "#E78AC3", putida = "#8DA0CB", ptotopetra = "#A6D854")
species_order <- names(species_cols)


# ------------------------------------------------------------------------------
# Samples, ordered by species
# ------------------------------------------------------------------------------

popmap  <- read.table(popmap_file, col.names = c("sample", "species"))
samples <- tibble(sample = read.table(fam_file)$V2) %>%
  mutate(fam_index = row_number(),
         species   = factor(popmap$species[match(sample, popmap$sample)], levels = species_order)) %>%
  arrange(species, sample) %>%
  mutate(order = row_number())


# ------------------------------------------------------------------------------
# Best replicate per K
# ------------------------------------------------------------------------------

last_number <- function(x) as.numeric(str_extract(x, "-?[0-9.]+$"))

admixture_best <- function(k) {
  cv <- map_dbl(seq_len(n_reps), function(r)
    last_number(read_lines(file.path(admix_dir, paste0("K", k), paste0("cv_K", k, "_rep", r, ".error")))[1]))
  best <- which.min(cv)
  message("ADMIXTURE K = ", k, ": rep ", best, ", CV = ", cv[best])
  as.matrix(read.table(file.path(admix_dir, paste0("K", k), paste0("K", k, "_rep", best, ".Q"))))
}

structure_q <- function(file, k) {
  lines <- read_lines(file)
  start <- grep("Inferred ancestry of individuals", lines)[1] + 2
  rows  <- lines[start:length(lines)]
  rows  <- rows[seq_len(which(str_trim(rows) == "")[1] - 1)]
  q     <- str_split(str_trim(str_remove(rows, ".*:")), "\\s+")
  do.call(rbind, map(q, ~ as.numeric(.x[seq_len(k)])))
}

structure_best <- function(k) {
  files <- file.path(struct_dir, paste0("K", k), paste0("K", k, "_rep", seq_len(n_reps), "_f"))
  lnp   <- map_dbl(files, ~ last_number(grep("Estimated Ln Prob of Data", read_lines(.x), value = TRUE)[1]))
  best  <- which.max(lnp)
  message("STRUCTURE K = ", k, ": rep ", best, ", lnP(D) = ", lnp[best])
  structure_q(files[best], k)
}


# ------------------------------------------------------------------------------
# Colours: each cluster takes the colour of the species with the highest mean
# membership; further clusters assigned to the same species are grey
# ------------------------------------------------------------------------------

cluster_colours <- function(q) {
  means <- apply(q, 2, function(x) tapply(x, samples$species[order(samples$fam_index)], mean))
  major <- species_order[apply(means, 2, which.max)]
  cols  <- species_cols[major]
  dup   <- duplicated(major[order(-apply(means, 2, max))])
  cols[order(-apply(means, 2, max))][dup] <- grey.colors(sum(dup), start = 0.3, end = 0.7)
  setNames(unname(cols), paste0("C", seq_len(ncol(q))))
}


# ------------------------------------------------------------------------------
# Bar plot for one K
# ------------------------------------------------------------------------------

bar_plot <- function(q, k, show_species = FALSE) {
  colnames(q) <- paste0("C", seq_len(k))
  cols <- cluster_colours(q)
  d <- as_tibble(q[samples$fam_index, , drop = FALSE]) %>%
    mutate(sample = samples$sample, order = samples$order) %>%
    pivot_longer(starts_with("C"), names_to = "cluster", values_to = "q")
  bounds <- samples %>% group_by(species) %>% summarise(end = max(order), mid = mean(order))

  p <- ggplot(d, aes(order, q, fill = cluster)) +
    geom_col(width = 1) +
    geom_vline(xintercept = head(bounds$end, -1) + 0.5, colour = "white", linewidth = 0.5) +
    scale_fill_manual(values = cols) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0), breaks = c(0, 0.5, 1)) +
    labs(x = NULL, y = paste0("K = ", k)) +
    theme_classic() +
    theme(legend.position = "none", axis.line.x = element_blank(),
          axis.ticks.x = element_blank(), axis.text.x = element_blank())
  if (show_species)
    p <- p + scale_x_continuous(expand = c(0, 0), breaks = bounds$mid,
                                labels = paste0("C. ", bounds$species)) +
      theme(axis.text.x = element_text(face = "italic", angle = 35, hjust = 1))
  p
}


# ------------------------------------------------------------------------------
# Fig. 1D (ADMIXTURE) and 1E (STRUCTURE)
# ------------------------------------------------------------------------------

admix_panels  <- map(k_values, ~ bar_plot(admixture_best(.x), .x, show_species = .x == max(k_values)))
struct_panels <- map(k_values, ~ bar_plot(structure_best(.x), .x, show_species = .x == max(k_values)))

fig <- wrap_plots(admix_panels, ncol = 1) | wrap_plots(struct_panels, ncol = 1)

ggsave(file.path(out_dir, "fig1_admixture_structure_K2-6.pdf"), fig, width = 11, height = 10)
