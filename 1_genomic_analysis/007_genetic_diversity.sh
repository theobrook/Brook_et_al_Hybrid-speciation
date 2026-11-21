#!/bin/bash
# ==============================================================================
# 007_genetic_diversity.sh
# Brook et al. (2026) Not a dead end: hybrid speciation and syngameon-driven
# dynamics in an island plant radiation
#
# Per-species observed heterozygosity (Ho), expected heterozygosity (He) and
# inbreeding coefficient (FIS), nucleotide diversity (pi), pairwise relative
# differentiation (FST) and absolute divergence (DXY)
# (Methods S2, "Genome-wide population genetic statistics"). Individuals with
# substantial mixed ancestry at K = 6 (C67, C70) were excluded.
#
# This file documents the commands used for each program. It is not a runnable
# pipeline: commands are shown for one example species and were run for each
# species.
#
# Software:
#   PLINK v.1.90b6.21
#   BCFtools v.1.21
#   VCFtools v.0.1.16
#   pixy v.2.0.0.beta14
# ==============================================================================


# ------------------------------------------------------------------------------
# Variables
# ------------------------------------------------------------------------------

PLINK_IN=snps                        # biallelic SNPs, no outgroups (006_population_structure.sh)
SPECIES=ptotopetra                   # example species
KEEPFILE=keep_files/${SPECIES}.keep  # two columns (FID IID) listing samples of this species
POPS_FILE=populations.txt            # two columns (sample, species), excluding C67 and C70
POP1=huttoniana                      # example species pair
POP2=ptotopetra


# ------------------------------------------------------------------------------
# 1. Remove individuals with mixed ancestry
# ------------------------------------------------------------------------------

printf "C67 C67\nC70 C70\n" > exclude_samples.txt
plink --bfile ${PLINK_IN} --remove exclude_samples.txt --make-bed --out snps_filtered


# ------------------------------------------------------------------------------
# 2. Heterozygosity and inbreeding coefficient per species
# ------------------------------------------------------------------------------

plink --bfile snps_filtered --keep ${KEEPFILE} --het   --out ${SPECIES}/het
plink --bfile snps_filtered --keep ${KEEPFILE} --hardy --out ${SPECIES}/hardy


# ------------------------------------------------------------------------------
# 3. All-sites VCF for pixy: biallelic SNPs (one per contig) + invariant sites
# ------------------------------------------------------------------------------
# Inputs from 005_vcf_filtering.sh. Also used for per-site FST outliers
# (009_selection_outliers.sh)

bcftools view --types snps --min-alleles 2 --max-alleles 2 \
    -Oz -o variant_sites_biallelic.vcf.gz variant_sites_stringent.vcf.gz
tabix -p vcf variant_sites_biallelic.vcf.gz

bcftools concat --allow-overlaps variant_sites_biallelic.vcf.gz invariant_sites.vcf.gz | \
    bcftools sort -Oz -o allsites_stringent.vcf.gz
tabix -p vcf allsites_stringent.vcf.gz


# ------------------------------------------------------------------------------
# 4. pixy: DXY for each species pair (run for all 28 pairs)
# ------------------------------------------------------------------------------
# Per-site output (window size 1); genome-wide DXY per pair was then calculated
# as sum(count_diffs) / sum(count_comparisons) (007_pixy_summary.R)

awk -v p1=${POP1} -v p2=${POP2} '$2 == p1 || $2 == p2' ${POPS_FILE} > ${POP1}_${POP2}_pops.txt

pixy --stats dxy \
    --vcf allsites_stringent.vcf.gz \
    --populations ${POP1}_${POP2}_pops.txt \
    --window_size 1 \
    --output_folder dxy \
    --output_prefix ${POP1}_vs_${POP2} \
    --n_cores 4


# ------------------------------------------------------------------------------
# 5. pixy: nucleotide diversity (pi) for each species
# ------------------------------------------------------------------------------
# Per-site output (window size 1); genome-wide pi per species was then
# calculated as sum(count_diffs) / sum(count_comparisons) (007_pixy_summary.R)

awk -v p=${SPECIES} '$2 == p' ${POPS_FILE} > ${SPECIES}_pops.txt

pixy --stats pi \
    --vcf allsites_stringent.vcf.gz \
    --populations ${SPECIES}_pops.txt \
    --window_size 1 \
    --output_folder pi \
    --output_prefix ${SPECIES}_pi \
    --n_cores 4


# ------------------------------------------------------------------------------
# 6. VCFtools: pairwise Weir & Cockerham FST between species
# ------------------------------------------------------------------------------
# Input: linkage-filtered variant sites (005_vcf_filtering.sh). One file per
# species listing its samples present in the VCF, excluding C67 and C70
# (pop_lists/<species>.txt). Run for every species pair; the weighted FST estimate (printed by VCFtools) is reported in Table S3
# and used in 007_fst_hybrid_test.R.

POP2_SPECIES=huttoniana               # example second species

vcftools --gzvcf variant_sites_stringent.vcf.gz \
    --weir-fst-pop pop_lists/${SPECIES}.txt \
    --weir-fst-pop pop_lists/${POP2_SPECIES}.txt \
    --out fst/${SPECIES}_vs_${POP2_SPECIES} \
    2> fst/${SPECIES}_vs_${POP2_SPECIES}.log

grep "Weir and Cockerham weighted Fst estimate:" fst/${SPECIES}_vs_${POP2_SPECIES}.log
