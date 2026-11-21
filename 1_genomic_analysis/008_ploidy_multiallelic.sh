#!/bin/bash
# ==============================================================================
# 008_ploidy_multiallelic.sh
# Brook et al. (2026) Not a dead end: hybrid speciation and syngameon-driven
# dynamics in an island plant radiation
#
# Multiallelic sites per individual, as a check for polyploidy (Methods S6;
# Table S16). Counted in the merged BCF before any filtering or linkage
# subsetting (004_variant_calling.sh). Analysed in 008_ploidy_multiallelic.R.
#
#   Total sites:        sites with a non-missing genotype in the individual
#   Multiallelic sites: of these, sites with >= 3 alleles at the locus
#
# This file documents the commands used. It is not a runnable pipeline.
#
# Software:
#   BCFtools v.1.21
# ==============================================================================


# ------------------------------------------------------------------------------
# Variables
# ------------------------------------------------------------------------------

MERGED_BCF=003_merged_output.bcf     # merged, unfiltered BCF (004_variant_calling.sh)
OUT_DIR=results/008_multiallelic


# ------------------------------------------------------------------------------
# 1. Sites with a non-missing genotype, per individual
# ------------------------------------------------------------------------------

bcftools query -f '[%SAMPLE\t%GT\n]' ${MERGED_BCF} \
    | awk '$2 != "./." && $2 != "." {n[$1]++} END {for (s in n) print s, n[s]}' \
    > ${OUT_DIR}/total_sites_per_sample.txt


# ------------------------------------------------------------------------------
# 2. Multiallelic sites (>= 3 alleles) with a non-missing genotype, per individual
# ------------------------------------------------------------------------------

bcftools view -m3 ${MERGED_BCF} \
    | bcftools query -f '[%SAMPLE\t%GT\n]' \
    | awk '$2 != "./." && $2 != "." {n[$1]++} END {for (s in n) print s, n[s]}' \
    > ${OUT_DIR}/multiallelic_sites_per_sample.txt


# ------------------------------------------------------------------------------
# 3. Percentage of multiallelic sites per individual (Table S16)
# ------------------------------------------------------------------------------

awk 'NR == FNR {tot[$1] = $2; next}
     {printf "%s\t%d\t%d\t%.4f\n", $1, tot[$1], $2, 100 * $2 / tot[$1]}' \
    ${OUT_DIR}/total_sites_per_sample.txt ${OUT_DIR}/multiallelic_sites_per_sample.txt \
    > ${OUT_DIR}/multiallelic_per_sample.tsv
