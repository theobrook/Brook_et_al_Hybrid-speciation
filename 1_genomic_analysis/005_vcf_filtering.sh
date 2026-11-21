#!/bin/bash
# ==============================================================================
# 005_vcf_filtering.sh
# Brook et al. (2026) Not a dead end: hybrid speciation and syngameon-driven
# dynamics in an island plant radiation
#
# Genotype and site filtering (Methods S1), and the linkage filter retaining one
# SNP per contig (main Methods, "Bioinformatic processing").
#
# This file documents the commands used for each program. It is not a runnable
# pipeline.
#
# Software:
#   BCFtools v.1.21
#   VCFtools v.0.1.16
# ==============================================================================


# ------------------------------------------------------------------------------
# Variables
# ------------------------------------------------------------------------------

BCF_IN=merged_tags.MaxDPfilt.bcf


# ------------------------------------------------------------------------------
# 1. Set genotypes with depth < 2 or genotype quality < 20 to missing
# ------------------------------------------------------------------------------

bcftools filter -S . -e 'FMT/GQ < 20 || FMT/DP < 2' ${BCF_IN} | \
    bgzip -c > basic_filt.vcf.gz
tabix -p vcf basic_filt.vcf.gz


# ------------------------------------------------------------------------------
# 2. Invariant sites (--max-maf 0 keeps only invariant sites)
# ------------------------------------------------------------------------------

vcftools --gzvcf basic_filt.vcf.gz \
    --max-maf 0 \
    --max-missing 0.80 \
    --minQ 20 \
    --min-meanDP 3 --max-meanDP 500 \
    --maxDP 500 \
    --recode --stdout | bgzip -c > invariant_sites.vcf.gz
tabix -p vcf invariant_sites.vcf.gz


# ------------------------------------------------------------------------------
# 3. Variant sites (--mac 1 keeps only variant sites; indels removed)
# ------------------------------------------------------------------------------

vcftools --gzvcf basic_filt.vcf.gz \
    --mac 1 \
    --remove-indels \
    --max-missing 0.80 \
    --minQ 20 \
    --min-meanDP 5 --max-meanDP 500 \
    --maxDP 500 \
    --recode --stdout | bgzip -c > variant_sites.vcf.gz
tabix -p vcf variant_sites.vcf.gz


# ------------------------------------------------------------------------------
# 4. Concatenate invariant and variant sites
# ------------------------------------------------------------------------------
# Combined all-sites dataset (no linkage filter)

bcftools concat -a -Ov invariant_sites.vcf.gz variant_sites.vcf.gz | bgzip -c > concat_sites.vcf.gz
bcftools index -t concat_sites.vcf.gz


# ------------------------------------------------------------------------------
# 5. Linkage filter: retain one SNP per contig (the site with lowest missingness)
# ------------------------------------------------------------------------------
# Used for population structure, genetic diversity and private alleles
# (n = 7,306 SNPs)

# Per-site missingness
vcftools --gzvcf variant_sites.vcf.gz --missing-site --out variant_sites_miss

# Select the SNP with the lowest missingness on each contig
tail -n +2 variant_sites_miss.lmiss | \
    sort -k1,1 -k6,6g | \
    awk '!seen[$1]++ {print $1"\t"$2}' > best_snp_per_contig.txt

# Keep only the selected SNPs
bcftools view -T best_snp_per_contig.txt variant_sites.vcf.gz -Oz -o variant_sites_stringent.vcf.gz
tabix -p vcf variant_sites_stringent.vcf.gz

# Combine with invariant sites
bcftools concat -a -Ov invariant_sites.vcf.gz variant_sites_stringent.vcf.gz | \
    bgzip -c > concat_sites_stringent.vcf.gz
bcftools index -t concat_sites_stringent.vcf.gz
