################################################################################
# 005_amf_composition.R
# Brook et al. (2026) Not a dead end: hybrid speciation and syngameon-driven
# dynamics in an island plant radiation
#
# Root AM fungal community composition among host species (Methods S4,
# "Analysis of AM fungal community composition").
#
#   1. Global PERMANOVA ~ host species (adonis2; Euclidean distance; 9,999
#      permutations) on Hellinger-transformed abundances and on presence/
#      absence, for ASVs and VTs, in batch 1 and in batches 1 + 2 (Table S11)
#   2. PCA of centred, Hellinger-transformed (or presence/absence) ASV and VT
#      matrices, batches 1 + 2 (Fig. S7)
#   3. Genus-level composition of soils and roots (Fig. S6)
#
# Inputs from 004_amf_tables.R. Only roots with >= 200 AM fungal reads are
# analysed. Pairwise and hybrid-parent contrasts were not performed (one to
# four roots per species).
#
# Software: R v.4.x, vegan v.2.7.1
################################################################################

library(tidyverse)
library(vegan)
library(patchwork)


# ------------------------------------------------------------------------------
# Settings
# ------------------------------------------------------------------------------

in_dir  <- "results/004_amf_tables"
out_dir <- "results/005_amf_composition"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

n_perm <- 9999
set.seed(42)

species_colours <- c(huttoniana = "#66C2A5", inopinata = "#FFD92F", lanceolaris = "#FC8D62",
                     ptotopetra = "#A6D854", putida = "#8DA0CB", savolainenii = "#E78AC3")

datasets <- c(batch1 = "Batch 1", batch1_2 = "Batches 1 + 2")


# ------------------------------------------------------------------------------
# Data: sample x taxon matrices
# ------------------------------------------------------------------------------

to_matrix <- function(d, taxon_col) {
  wide <- unclass(xtabs(as.formula(paste("reads ~ Sample_ID +", taxon_col)), data = d))
  meta <- d %>% distinct(Sample_ID, Species, Batch) %>%
    mutate(Species = factor(Species), Batch = factor(Batch))
  meta <- meta[match(rownames(wide), meta$Sample_ID), ]
  list(abundance = decostand(wide, "hellinger"),
       presence_absence = (wide > 0) * 1,
       meta = meta)
}

matrices <- list()
for (ds in names(datasets)) {
  matrices[[paste("ASV", ds)]] <- to_matrix(read_csv(file.path(in_dir, paste0("asv_long_amf_", ds, ".csv")),
                                                     show_col_types = FALSE), "asv")
  matrices[[paste("VT", ds)]]  <- to_matrix(read_csv(file.path(in_dir, paste0("vt_long_amf_", ds, ".csv")),
                                                     show_col_types = FALSE), "vt_code")
}


# ------------------------------------------------------------------------------
# 1. Global PERMANOVA (Table S11)
# ------------------------------------------------------------------------------

permanova <- imap_dfr(matrices, function(x, nm) {
  map_dfr(c("abundance", "presence_absence"), function(tr) {
    a <- adonis2(x[[tr]] ~ Species, data = x$meta, method = "euclidean",
                 permutations = n_perm)
    tibble(level     = word(nm, 1),
           dataset   = datasets[[word(nm, 2)]],
           data      = tr,
           n_roots   = nrow(x[[tr]]),
           n_species = nlevels(x$meta$Species),
           df        = paste(a$Df[1], a$Df[2], sep = ", "),
           F         = a$F[1],
           R2        = a$R2[1],
           p         = a$`Pr(>F)`[1])
  })
})

write_csv(permanova, file.path(out_dir, "global_permanova.csv"))


# ------------------------------------------------------------------------------
# 2. PCA, batches 1 + 2 (Fig. S7)
# ------------------------------------------------------------------------------

pca_panel <- function(m, meta, title) {
  pca    <- prcomp(m, center = TRUE, scale. = FALSE)
  var    <- round(100 * summary(pca)$importance[2, 1:2], 1)
  scores <- as_tibble(pca$x[, 1:2]) %>% bind_cols(meta)

  ggplot(scores, aes(PC1, PC2, fill = Species, shape = Batch)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey70", linewidth = 0.3) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey70", linewidth = 0.3) +
    geom_point(size = 3, colour = "white", stroke = 0.3) +
    scale_fill_manual(values = species_colours,
                      labels = function(x) parse(text = paste0("italic('C. ", x, "')"))) +
    scale_shape_manual(values = c("1" = 21, "2" = 24)) +
    guides(fill = guide_legend(override.aes = list(shape = 21))) +
    labs(title = title, x = paste0("PC1 (", var[1], "%)"), y = paste0("PC2 (", var[2], "%)")) +
    theme_minimal(base_size = 11) +
    theme(panel.grid.minor = element_blank(),
          panel.border = element_rect(colour = "grey80", fill = NA))
}

asv <- matrices[["ASV batch1_2"]]
vt  <- matrices[["VT batch1_2"]]

fig_s7 <- pca_panel(asv$abundance,        asv$meta, "A  ASV abundance") +
          pca_panel(asv$presence_absence, asv$meta, "B  ASV presence/absence") +
          pca_panel(vt$abundance,         vt$meta,  "C  VT abundance") +
          pca_panel(vt$presence_absence,  vt$meta,  "D  VT presence/absence") +
          plot_layout(ncol = 2, guides = "collect")

ggsave(file.path(out_dir, "fig_s7_root_amf_pca.png"), fig_s7,
       width = 12, height = 10, dpi = 300, bg = "white")


# ------------------------------------------------------------------------------
# 3. Genus composition of soils and roots (Fig. S6)
# ------------------------------------------------------------------------------
# Non-rarefied AM fungal reads; bars labelled with sample number and AM fungal
# reads per sample; roots included in the analyses (>= 200 AM fungal reads)
# marked with *

amf_reads <- read_csv(file.path(in_dir, "amf_reads_per_sample.csv"), show_col_types = FALSE)
comp <- read_csv(file.path(in_dir, "asv_long_amf_nonrarefied.csv"), show_col_types = FALSE) %>%
  mutate(genus = replace_na(genus, "Unassigned")) %>%
  group_by(Sample_ID, Species, Sample_type, genus) %>%
  summarise(reads = sum(reads), .groups = "drop") %>%
  group_by(Sample_ID) %>%
  mutate(rel_abund = reads / sum(reads)) %>%
  ungroup()

composition_panel <- function(st, title) {
  samp <- amf_reads %>%
    filter(Sample_type == st) %>%
    mutate(pair  = str_remove(Sample_ID, "[RS]$"),
           label = paste0(pair, " (", amf_reads, ")", if_else(retained, "*", ""))) %>%
    arrange(Species, as.numeric(str_extract(pair, "[0-9]+")))
  d <- comp %>% filter(Sample_type == st) %>%
    mutate(Sample_ID = factor(Sample_ID, levels = intersect(samp$Sample_ID, Sample_ID)))

  ggplot(d, aes(Sample_ID, rel_abund, fill = genus)) +
    geom_col(width = 0.85) +
    facet_grid(~ Species, scales = "free_x", space = "free_x") +
    scale_x_discrete(labels = setNames(samp$label, samp$Sample_ID)) +
    scale_y_continuous(expand = c(0, 0)) +
    labs(title = title, x = NULL, y = "Relative abundance", fill = "Genus") +
    theme_classic(base_size = 10) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
          strip.text = element_text(face = "italic"),
          strip.background = element_blank())
}

fig_s6 <- composition_panel("soil", "A  Soil") /
          composition_panel("root", "B  Root") +
          plot_layout(guides = "collect")

ggsave(file.path(out_dir, "fig_s6_amf_genus_composition.png"), fig_s6,
       width = 12, height = 9, dpi = 300, bg = "white")
