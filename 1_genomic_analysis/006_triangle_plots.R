################################################################################
# 006_triangle_plots.R
# Brook et al. (2026) Not a dead end: hybrid speciation and syngameon-driven
# dynamics in an island plant radiation
#
# Triangle plots of hybrid index against interclass heterozygosity, used to
# distinguish stabilised hybrid lineages from early-generation hybrids
# (Methods S2, "Triangle plots"; Fig. 2).
#
# Input: biallelic, linkage-filtered SNPs without outgroups
#        (biallelic_snps.vcf.gz from 006_population_structure.sh) and a
#        population map (popmap.csv; columns: id, pop).
#
# Parental combinations tested:
#   A  C. putida        = C. ptotopetra x C. huttoniana
#   B  C. savolainenii  = C. ptotopetra x C. inopinata
#   C  C. savolainenii  = C. ptotopetra x C. lanceolaris
# Ancestry-informative markers (AIMs): allele frequency difference >= 0.9
# between the two parental species.
################################################################################

library(triangulaR)
library(vcfR)
library(tidyverse)
library(patchwork)


# ------------------------------------------------------------------------------
# Settings
# ------------------------------------------------------------------------------

vcf_file    <- "biallelic_snps.vcf.gz"
popmap_file <- "popmap.csv"
out_dir     <- "results/triangle_plots"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

aim_threshold <- 0.9
recent_hybrids <- c(C67 = 3, C70 = 8)   # highlighted with point shapes 3 (+) and 8 (*)

species_colours <- c(
  putida       = "#8DA0CB",
  lanceolaris  = "#FC8D62",
  ptotopetra   = "#A6D854",
  huttoniana   = "#66C2A5",
  inopinata    = "#FFD92F",
  savolainenii = "#E78AC3"
)


# ------------------------------------------------------------------------------
# Data
# ------------------------------------------------------------------------------

vcf    <- read.vcfR(vcf_file)
popmap <- read.csv(popmap_file)


# ------------------------------------------------------------------------------
# Functions
# ------------------------------------------------------------------------------

# Subset the VCF and population map to the two parents and the hybrid species
subset_species <- function(vcf, popmap, keep) {
  pm   <- popmap %>% filter(pop %in% keep)
  cols <- which(colnames(vcf@gt)[-1] %in% pm$id) + 1
  vcf@gt <- vcf@gt[, c(1, cols), drop = FALSE]
  list(vcf = vcf, popmap = pm)
}

# Hybrid index and interclass heterozygosity from AIMs, then plot
triangle_plot <- function(vcf, popmap, p1, p2, hybrid) {

  d      <- subset_species(vcf, popmap, c(p1, p2, hybrid))
  aims   <- alleleFreqDiff(vcfR = d$vcf, pm = d$popmap, p1 = p1, p2 = p2,
                           difference = aim_threshold)
  hi_het <- hybridIndex(vcfR = aims, pm = d$popmap, p1 = p1, p2 = p2)

  hwe <- tibble(x = seq(0, 1, length.out = 500), y = 2 * x * (1 - x))
  rec <- hi_het %>% filter(id %in% names(recent_hybrids))

  ggplot(hi_het, aes(hybrid.index, heterozygosity)) +
    annotate("segment", x = 0,   xend = 0.5, y = 0, yend = 1) +
    annotate("segment", x = 0.5, xend = 1,   y = 1, yend = 0) +
    annotate("segment", x = 0,   xend = 1,   y = 0, yend = 0) +
    geom_path(data = hwe, aes(x, y), linetype = "dashed", inherit.aes = FALSE) +
    geom_point(aes(fill = pop), shape = 21, size = 6, colour = "transparent") +
    geom_point(data = rec, aes(shape = id), colour = species_colours[["putida"]],
               size = 7, stroke = 0.7) +
    scale_fill_manual(values = species_colours) +
    scale_shape_manual(values = recent_hybrids) +
    coord_cartesian(xlim = c(-0.02, 1.02), ylim = c(-0.02, 1.02)) +
    labs(x = "Hybrid Index (HI)", y = "Interclass heterozygosity") +
    theme_classic() +
    theme(legend.position = "none",
          axis.text  = element_text(size = 15, colour = "black"),
          axis.title = element_text(size = 18, colour = "black"),
          plot.tag   = element_text(size = 18, face = "bold"))
}


# ------------------------------------------------------------------------------
# Triangle plots
# ------------------------------------------------------------------------------

p_putida      <- triangle_plot(vcf, popmap, "ptotopetra", "huttoniana",  "putida")
p_savo_inop   <- triangle_plot(vcf, popmap, "ptotopetra", "inopinata",   "savolainenii")
p_savo_lanc   <- triangle_plot(vcf, popmap, "ptotopetra", "lanceolaris", "savolainenii")

combined <- (p_putida | p_savo_inop | p_savo_lanc) +
  plot_annotation(tag_levels = "A")

ggsave(file.path(out_dir, "triangle_plots.pdf"), combined, width = 18, height = 6.5)
