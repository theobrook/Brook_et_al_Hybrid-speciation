################################################################################
# 004_amf_tables.R
# Brook et al. (2026) Not a dead end: hybrid speciation and syngameon-driven
# dynamics in an island plant radiation
#
# AM fungal ASV and virtual taxon (VT) tables from the BLAST searches in
# 003_amf_identification.sh (Methods S4, "Analysis of AM fungal community
# composition"; Tables S11 and S19).
#
#   1. MaarjAM: best hit per ASV by bitscore; VT code, family and genus parsed
#      from the hit title
#   2. SILVA check: an ASV is kept as AM fungal only if all of its best-scoring
#      SILVA hits are Glomeromycota
#   3. AM fungal reads per sample, taken from the non-rarefied ASV table
#   4. Roots with < 200 AM fungal reads excluded; AM fungal reads of the
#      retained roots rarefied (rrarefy) to the smallest depth among them;
#      VT tables built by summing the reads of ASVs assigned to the same VT.
#      Done for batch 1 and for batches 1 + 2.
#
# Software: R v.4.x, vegan v.2.7.1
################################################################################

library(tidyverse)
library(vegan)


# ------------------------------------------------------------------------------
# Settings
# ------------------------------------------------------------------------------

blast_dir     <- "results/003_blast"
dada2_dir     <- "results/002_dada2"
out_dir       <- "results/004_amf_tables"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

min_amf_reads <- 200      # roots: minimum AM fungal reads
rare_seed     <- 1

blast_cols <- c("asv", "hit", "pident", "length", "qlen", "qcovs", "evalue", "bitscore", "stitle")
glom_pattern <- paste("Glomeromycot", "Glomeromycetes", "Glomerales", "Paraglomerales",
                      "Archaeosporales", "Diversisporales", "Glomeraceae", "Paraglomeraceae",
                      sep = "|")


# ------------------------------------------------------------------------------
# 1. MaarjAM best hit per ASV
# ------------------------------------------------------------------------------

maarjam <- read_tsv(file.path(blast_dir, "blast_maarjam.txt"), col_names = blast_cols,
                    show_col_types = FALSE) %>%
  group_by(asv) %>%
  slice_max(bitscore, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(vt_code = str_replace(stitle, ".*(VTX[0-9]+).*", "\\1"),
         family  = str_replace(stitle, "^gi\\|[^|]+\\|gb\\|[^|]+\\| ([A-Za-z]+).*", "\\1"),
         genus   = str_replace(stitle, "^gi\\|[^|]+\\|gb\\|[^|]+\\| [A-Za-z]+ ([A-Za-z]+).*", "\\1")) %>%
  filter(str_detect(vt_code, "^VTX[0-9]+$"))


# ------------------------------------------------------------------------------
# 2. SILVA check
# ------------------------------------------------------------------------------

silva <- read_tsv(file.path(blast_dir, "blast_silva.txt"), col_names = blast_cols,
                  show_col_types = FALSE, quote = "") %>%
  mutate(taxonomy = str_remove(stitle, "^\\S+\\s+"),
         is_glom  = str_detect(taxonomy, glom_pattern))

silva_verdict <- silva %>%
  group_by(asv) %>%
  filter(bitscore == max(bitscore)) %>%
  summarise(all_glom = all(is_glom), .groups = "drop")

amf_asvs <- maarjam %>%
  inner_join(silva_verdict, by = "asv") %>%
  filter(all_glom) %>%
  select(asv, vt_code, family, genus, pident, qcovs)

write_csv(amf_asvs, file.path(out_dir, "amf_asvs.csv"))


# ------------------------------------------------------------------------------
# 3. AM fungal reads per sample (non-rarefied ASV table)
# ------------------------------------------------------------------------------

nonrare <- read_csv(file.path(dada2_dir, "asv_long_nonrarefied.csv"), show_col_types = FALSE) %>%
  filter(Species != "prisca")

amf_reads <- nonrare %>%
  group_by(Sample_ID, Species, Sample_type, Batch) %>%
  summarise(total_reads = sum(reads),
            amf_reads   = sum(reads[asv %in% amf_asvs$asv]), .groups = "drop") %>%
  mutate(retained = Sample_type == "root" & amf_reads >= min_amf_reads)

write_csv(amf_reads, file.path(out_dir, "amf_reads_per_sample.csv"))

# Non-rarefied AM fungal reads of all samples (genus composition, Fig. S6)
nonrare %>%
  filter(asv %in% amf_asvs$asv) %>%
  left_join(amf_asvs, by = "asv") %>%
  write_csv(file.path(out_dir, "asv_long_amf_nonrarefied.csv"))


# ------------------------------------------------------------------------------
# 4. Retained roots, rarefied AM fungal ASV and VT tables, per dataset
# ------------------------------------------------------------------------------

datasets <- list(batch1 = 1, batch1_2 = c(1, 2))

for (ds in names(datasets)) {

  roots <- amf_reads %>% filter(retained, Batch %in% datasets[[ds]])

  amf_wide <- nonrare %>%
    filter(asv %in% amf_asvs$asv, Sample_ID %in% roots$Sample_ID) %>%
    xtabs(reads ~ Sample_ID + asv, data = .) %>%
    unclass()

  set.seed(rare_seed)
  rare <- rrarefy(amf_wide, sample = min(rowSums(amf_wide)))
  rare <- rare[, colSums(rare) > 0, drop = FALSE]

  i <- which(rare > 0, arr.ind = TRUE)
  asv_long <- tibble(Sample_ID = rownames(rare)[i[, 1]],
                     asv       = colnames(rare)[i[, 2]],
                     reads     = rare[i]) %>%
    left_join(roots %>% select(Sample_ID, Species, Sample_type, Batch), by = "Sample_ID") %>%
    left_join(amf_asvs, by = "asv")

  vt_long <- asv_long %>%
    group_by(Sample_ID, Species, Sample_type, Batch, vt_code) %>%
    summarise(reads = sum(reads), family = first(family), genus = first(genus),
              .groups = "drop")

  write_csv(asv_long, file.path(out_dir, paste0("asv_long_amf_", ds, ".csv")))
  write_csv(vt_long,  file.path(out_dir, paste0("vt_long_amf_",  ds, ".csv")))

  message(ds, ": ", nrow(rare), " roots rarefied to ", min(rowSums(amf_wide)),
          " AM fungal reads; ", ncol(rare), " ASVs; ", n_distinct(vt_long$vt_code), " VTs")
}
