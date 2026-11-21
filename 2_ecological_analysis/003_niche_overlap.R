################################################################################
# 003_niche_overlap.R
# Brook et al. (2026) Not a dead end: hybrid speciation and syngameon-driven
# dynamics in an island plant radiation
#
# Pairwise niche overlap from the species distribution models (Methods S3,
# "Species distribution analysis"; Table S9). Run after 003_maxent.sh.
#
#   1. MaxEnt habitat suitability surface for each species, averaged across
#      cross-validation replicates (<species>_avg.asc)
#   2. Each surface normalised so that cell values sum to one
#   3. Schoener's D and Warren's I for all 15 species pairs:
#        D = 1 - 0.5 * sum(|p1 - p2|)
#        I = 1 - 0.5 * sqrt(sum((sqrt(p1) - sqrt(p2))^2))
#      Both range from 0 (no overlap) to 1 (identical distributions)
#
# The matrices are used in 004_fst_niche_correlations.R.
#
# Software: R v.4.x, raster v.3.6.32
################################################################################

library(raster)


# ------------------------------------------------------------------------------
# Settings
# ------------------------------------------------------------------------------

maxent_dirs <- c("maxent/no_thinning", "maxent/no_thinning_inopinata")  # 003_maxent.sh
out_dir     <- "results"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

focal_species <- c("huttoniana", "inopinata", "lanceolaris",
                   "ptotopetra", "putida", "savolainenii")


# ------------------------------------------------------------------------------
# 1. Averaged suitability surfaces
# ------------------------------------------------------------------------------

suitability <- lapply(focal_species, function(sp) {
  f <- file.path(maxent_dirs, paste0(sp, "_avg.asc"))
  raster(f[file.exists(f)][1])
})
names(suitability) <- focal_species


# ------------------------------------------------------------------------------
# 2-3. Normalised surfaces, Schoener's D and Warren's I
# ------------------------------------------------------------------------------

normalise <- function(p) p / sum(p)

niche_overlap <- function(r1, r2) {
  p1 <- values(r1)
  p2 <- values(r2)
  valid <- !is.na(p1) & !is.na(p2)
  p1 <- normalise(p1[valid])
  p2 <- normalise(p2[valid])
  c(D = 1 - 0.5 * sum(abs(p1 - p2)),
    I = 1 - 0.5 * sqrt(sum((sqrt(p1) - sqrt(p2))^2)))
}

n <- length(focal_species)
D_matrix <- I_matrix <- matrix(1, n, n, dimnames = list(focal_species, focal_species))

for (i in seq_len(n - 1)) {
  for (j in (i + 1):n) {
    ov <- niche_overlap(suitability[[i]], suitability[[j]])
    D_matrix[i, j] <- D_matrix[j, i] <- ov[["D"]]
    I_matrix[i, j] <- I_matrix[j, i] <- ov[["I"]]
  }
}

write.csv(round(D_matrix, 4), file.path(out_dir, "schoeners_D_matrix.csv"))
write.csv(round(I_matrix, 4), file.path(out_dir, "warrens_I_matrix.csv"))
