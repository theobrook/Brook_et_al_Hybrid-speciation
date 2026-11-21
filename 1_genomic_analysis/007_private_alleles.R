################################################################################
# 007_private_alleles.R
# Brook et al. (2026) Not a dead end: hybrid speciation and syngameon-driven
# dynamics in an island plant radiation
#
# Private alleles and fixed private differences per species, with
# individual-based rarefaction (Methods S2, "Private alleles and fixed private
# alleles").
#
# Input: biallelic, linkage-filtered SNPs without outgroups
#        (biallelic_snps.vcf.gz from 006_population_structure.sh) and a
#        two-column population map (sample, species).
#
#   1. Exclude individuals with mixed ancestry at K = 6 (C67, C70)
#   2. Per species: ALT allele frequency, proportion missing, number genotyped
#   3. Classify private alleles and fixed private differences:
#        focal species:      <= 20% missing genotypes
#        non-focal species:  >= 3 genotyped individuals each
#        private:            allele present in focal species, absent elsewhere
#        fixed private:      allele fixed in focal species, absent elsewhere
#      The REF allele comes from the assembly sample (C11, C. lanceolaris), so
#      it is known to occur in C. lanceolaris. REF can therefore only be private
#      (or fixed private) in C. lanceolaris, and C. lanceolaris cannot be fixed
#      for ALT.
#   4. Rarefaction: 1,000 iterations subsampling each species to the smallest
#      sample size (n = 5) and recounting
################################################################################

library(vcfR)
library(tidyverse)


# ------------------------------------------------------------------------------
# Settings
# ------------------------------------------------------------------------------

vcf_file    <- "biallelic_snps.vcf.gz"
popmap_file <- "samples.txt"
output_dir  <- "results/private_alleles"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

max_missing_focal <- 0.20
min_nn_other      <- 3
mixed_samples     <- c("C67", "C70")
reference_species <- "lanceolaris"
n_iterations      <- 1000


# ------------------------------------------------------------------------------
# 1. Load genotypes as ALT allele counts (0/1/2, NA = missing)
# ------------------------------------------------------------------------------

vcf    <- read.vcfR(vcf_file, verbose = FALSE)
popmap <- read.table(popmap_file, col.names = c("sample_id", "species")) %>%
  filter(sample_id %in% colnames(vcf@gt)[-1], !sample_id %in% mixed_samples)

species_list <- unique(popmap$species)

gt <- extract.gt(vcf, element = "GT")
gt <- apply(gt, 2, function(x) {
  sapply(strsplit(gsub("\\|", "/", x), "/"), function(a) {
    if (any(is.na(a)) || any(a == ".")) NA else sum(as.numeric(a))
  })
})
gt <- gt[, popmap$sample_id, drop = FALSE]


# ------------------------------------------------------------------------------
# 2. Per-species allele frequency, missingness and number genotyped
# ------------------------------------------------------------------------------

species_stats <- function(gt, pm) {
  out <- lapply(unique(pm$species), function(sp) {
    g  <- gt[, pm$sample_id[pm$species == sp], drop = FALSE]
    nn <- rowSums(!is.na(g))
    list(af   = ifelse(nn > 0, rowSums(g, na.rm = TRUE) / (2 * nn), NA),
         miss = 1 - nn / ncol(g),
         nn   = nn)
  })
  setNames(out, unique(pm$species))
}


# ------------------------------------------------------------------------------
# 3. Classify private alleles (fixed = FALSE) or fixed private differences
#    (fixed = TRUE) for a focal species
# ------------------------------------------------------------------------------

classify_private <- function(st, focal, fixed = FALSE) {
  others   <- setdiff(names(st), focal)
  f        <- st[[focal]]
  focal_ok <- !is.na(f$af) & f$miss <= max_missing_focal

  others_ok  <- Reduce(`&`, lapply(others, function(sp) !is.na(st[[sp]]$af) & st[[sp]]$nn >= min_nn_other))
  others_ref <- Reduce(`&`, lapply(others, function(sp) st[[sp]]$af == 0))
  others_alt <- Reduce(`&`, lapply(others, function(sp) st[[sp]]$af == 1))

  is_ref_sp <- focal == reference_species
  none      <- rep(FALSE, length(f$af))

  if (fixed) {
    alt <- if (is_ref_sp) none else focal_ok & f$af == 1 & others_ok & others_ref
    ref <- if (is_ref_sp) focal_ok & f$af == 0 & others_ok & others_alt else none
  } else {
    alt <- focal_ok & f$af > 0 & others_ok & others_ref
    ref <- if (is_ref_sp) focal_ok & f$af < 1 & others_ok & others_alt else none
  }
  (alt | ref) %in% TRUE
}

st <- species_stats(gt, popmap)

private_summary <- tibble(
  Species           = species_list,
  N_samples         = map_int(species_list, ~ sum(popmap$species == .x)),
  N_private_alleles = map_int(species_list, ~ sum(classify_private(st, .x, fixed = FALSE))),
  N_fixed_private   = map_int(species_list, ~ sum(classify_private(st, .x, fixed = TRUE)))
) %>%
  mutate(Prop_private_pct = 100 * N_private_alleles / nrow(gt),
         Prop_fixed_pct   = 100 * N_fixed_private   / nrow(gt))

write_tsv(private_summary, file.path(output_dir, "private_alleles_summary.tsv"))


# ------------------------------------------------------------------------------
# 4. Rarefaction to the smallest species sample size
# ------------------------------------------------------------------------------

set.seed(42)
rarefaction_n <- min(table(popmap$species))

rarefied <- map_dfr(seq_len(n_iterations), function(iter) {
  keep <- unlist(lapply(species_list, function(sp)
    sample(popmap$sample_id[popmap$species == sp], rarefaction_n)))
  sub    <- popmap %>% filter(sample_id %in% keep)
  st_sub <- species_stats(gt[, sub$sample_id, drop = FALSE], sub)
  tibble(iteration     = iter,
         Species       = species_list,
         private       = map_int(species_list, ~ sum(classify_private(st_sub, .x, FALSE))),
         fixed_private = map_int(species_list, ~ sum(classify_private(st_sub, .x, TRUE))))
})

rarefaction_summary <- rarefied %>%
  pivot_longer(c(private, fixed_private), names_to = "category", values_to = "count") %>%
  group_by(Species, category) %>%
  summarise(Rarefied_n  = rarefaction_n,
            Mean        = mean(count),
            Median      = median(count),
            SD          = sd(count),
            CI_lower_95 = quantile(count, 0.025),
            CI_upper_95 = quantile(count, 0.975),
            .groups = "drop")

write_tsv(rarefied,            file.path(output_dir, "rarefaction_iterations.tsv"))
write_tsv(rarefaction_summary, file.path(output_dir, "rarefaction_summary.tsv"))
