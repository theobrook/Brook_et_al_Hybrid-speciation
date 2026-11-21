################################################################################
# 009_fst_outliers.R
# Brook et al. (2026) Not a dead end: hybrid speciation and syngameon-driven
# dynamics in an island plant radiation
#
# Identifies per-site FST outliers between species pairs from pixy output and
# maps outlier sites to transcripts for annotation (main Methods, "Selection
# analyses"). Run after step 1 of 009_selection_outliers.sh.
#
#   1. Load per-site Weir-Cockerham FST for all species pairs (pixy, window = 1)
#   2. Remove outgroup comparisons (C. prisca, C. repens): 15 pairs remain
#   3. Flag the top 10% and top 5% of FST values in each pair. Quantile type 1
#      is used so that all sites tied at a threshold (e.g. FST = 1) are retained
#   4. Map outlier sites to transcript IDs and write the transcript lists used
#      by scripts/annotate_outliers_eggnog.py
################################################################################

library(tidyverse)
library(data.table)


# ------------------------------------------------------------------------------
# Paths
# ------------------------------------------------------------------------------

fst_dir     <- "results/outliers/fst"                          # pixy output
outlier_dir <- "results/outliers/002_fst_outliers"             # per-pair outlier tables
annot_dir   <- "results/outliers/003_annotation/pairwise_outliers"  # transcript lists

dir.create(outlier_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(annot_dir,   recursive = TRUE, showWarnings = FALSE)

outgroups <- c("prisca", "repens")


# ------------------------------------------------------------------------------
# 1. Load per-site FST
# ------------------------------------------------------------------------------

fst_files <- list.files(fst_dir, pattern = "_fst\\.txt$",
                        recursive = TRUE, full.names = TRUE)

fst <- rbindlist(lapply(fst_files, fread, na.strings = c("nan", "NA", ""))) %>%
  rename(Pop1 = pop1, Pop2 = pop2, Chrom = chromosome,
         Start = window_pos_1, End = window_pos_2) %>%
  filter(!Pop1 %in% outgroups, !Pop2 %in% outgroups) %>%
  filter(!is.na(avg_wc_fst)) %>%
  mutate(chrom_pos = paste(Chrom, Start, sep = "_")) %>%
  distinct()

species_pairs <- fst %>% distinct(Pop1, Pop2)


# ------------------------------------------------------------------------------
# 2. Flag outliers in each species pair
# ------------------------------------------------------------------------------

flag_outliers <- function(df) {
  q90 <- quantile(df$avg_wc_fst, 0.90, type = 1)
  q95 <- quantile(df$avg_wc_fst, 0.95, type = 1)
  df %>%
    arrange(Chrom, Start) %>%
    mutate(FST_outlier_top10 = avg_wc_fst >= q90,
           FST_outlier_top5  = avg_wc_fst >= q95,
           threshold_90      = q90,
           threshold_95      = q95)
}

fst_flagged <- fst %>%
  group_by(Pop1, Pop2) %>%
  group_modify(~ flag_outliers(.x)) %>%
  ungroup()


# ------------------------------------------------------------------------------
# 3. Summary across pairs
# ------------------------------------------------------------------------------

outlier_summary <- fst_flagged %>%
  group_by(Pop1, Pop2) %>%
  summarise(n_sites          = n(),
            n_sites_Fst1     = sum(avg_wc_fst == 1),
            n_sites_Fst0     = sum(avg_wc_fst == 0),
            FST_threshold_90 = first(threshold_90),
            FST_threshold_95 = first(threshold_95),
            n_outliers_top10 = sum(FST_outlier_top10),
            n_outliers_top5  = sum(FST_outlier_top5),
            .groups = "drop")

write_csv(outlier_summary, file.path(outlier_dir, "FST_outliers_summary.csv"))


# ------------------------------------------------------------------------------
# 4. Map outlier sites to transcripts and write lists for annotation
# ------------------------------------------------------------------------------
# Site IDs are "<contig>_<position>"; transcript (peptide) IDs are "<contig>.p1"

fst_flagged <- fst_flagged %>%
  mutate(transcript_id = paste0(sub("_[0-9]+$", "", chrom_pos), ".p1"))

for (i in seq_len(nrow(species_pairs))) {

  pop1 <- species_pairs$Pop1[i]
  pop2 <- species_pairs$Pop2[i]
  pair <- paste0(pop1, "_", pop2)

  pair_df <- fst_flagged %>% filter(Pop1 == pop1, Pop2 == pop2)
  write_csv(pair_df, file.path(outlier_dir, paste0(pair, "_FST_outliers.csv")))

  all_transcripts <- unique(pair_df$transcript_id)

  for (metric in c("FST_top10", "FST_top5")) {
    flag_col <- ifelse(metric == "FST_top10", "FST_outlier_top10", "FST_outlier_top5")
    out <- file.path(annot_dir, pair, metric)
    dir.create(out, recursive = TRUE, showWarnings = FALSE)

    write_lines(all_transcripts, file.path(out, "all_transcripts.txt"))
    write_lines(unique(pair_df$transcript_id[pair_df[[flag_col]]]),
                file.path(out, "outlier_transcripts.txt"))
  }
}
