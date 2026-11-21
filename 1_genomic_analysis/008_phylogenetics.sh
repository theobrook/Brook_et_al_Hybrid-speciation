#!/bin/bash
# ==============================================================================
# 008_phylogenetics.sh
# Brook et al. (2026) Not a dead end: hybrid speciation and syngameon-driven
# dynamics in an island plant radiation
#
# Individual-level maximum likelihood phylogeny (main Methods, "Bioinformatic
# processing"; Methods S2, "Individual-level phylogeny"). Includes three
# outgroups: C. prisca (C4, C43) and C. repens (C77).
#
# This file documents the commands used for each program. It is not a runnable
# pipeline: per-sample commands are shown for one example sample and were run
# for every sample.
#
# Software:
#   BCFtools v.1.21
#   VCFtools v.0.1.16
#   vcf2fas
#   stringent_filter_fasta_alns.py and cat_alns.py (scripts/)
#   RAxML-NG v.1.2.2
# ==============================================================================


# ------------------------------------------------------------------------------
# Variables
# ------------------------------------------------------------------------------

REFERENCE=Trinity.SuperTrans.fasta
BCF_DIR=results/011_bcftools_call      # per-sample BCFs (004_variant_calling.sh)

SAMPLE=C10                             # example sample

# All samples in the phylogeny (C11 was used only for the reference assembly)
TAXA="C10,C14,C16,C17,C18,C19,C21,C22,C24,C25,C26,C27,C2,C31,C32,C33,C34,C35,C36,C37,C41,C42,C45,C46,C48,C50,C53,C55,C56,C5,C8,C60,C67,C69,C70,C71,C72,C73,C74,C58,C79,C80,C81,C82,C83,C84,C4,C43,C77"


# ------------------------------------------------------------------------------
# 1. Per-sample GVCFs: same filters as for the multi-sample dataset
#    (005_vcf_filtering.sh), with no linkage filter
# ------------------------------------------------------------------------------

# Add tags and remove sites with excessive coverage
bcftools +fill-tags ${BCF_DIR}/${SAMPLE}call.bcf -Ob -o ${SAMPLE}_tags.bcf -- -t all
bcftools view -i 'INFO/DP < 1260' ${SAMPLE}_tags.bcf -Ob -o ${SAMPLE}_tags.MaxDPfilt.bcf

# Set genotypes with depth < 2 or genotype quality < 20 to missing
bcftools filter -S . -e 'FMT/GQ < 20 || FMT/DP < 2' ${SAMPLE}_tags.MaxDPfilt.bcf | \
    bgzip -c > ${SAMPLE}.basicfilt.vcf.gz
tabix -p vcf ${SAMPLE}.basicfilt.vcf.gz

# Invariant sites
vcftools --gzvcf ${SAMPLE}.basicfilt.vcf.gz \
    --max-maf 0 --max-missing 0.80 --minQ 20 \
    --min-meanDP 3 --max-meanDP 500 --maxDP 500 \
    --recode --stdout | bgzip -c > ${SAMPLE}.invariant.vcf.gz
tabix -p vcf ${SAMPLE}.invariant.vcf.gz

# Variant sites
vcftools --gzvcf ${SAMPLE}.basicfilt.vcf.gz \
    --mac 1 --remove-indels --max-missing 0.80 --minQ 20 \
    --min-meanDP 5 --max-meanDP 500 --maxDP 500 \
    --recode --stdout | bgzip -c > ${SAMPLE}.variant.vcf.gz
tabix -p vcf ${SAMPLE}.variant.vcf.gz

# Combine invariant and variant sites
bcftools concat -a -Ov ${SAMPLE}.invariant.vcf.gz ${SAMPLE}.variant.vcf.gz -o ${SAMPLE}.concat.vcf

# Number of variant sites per sample after filtering
grep -v "^#" ${SAMPLE}.concat.vcf | awk '$5 != "."' | wc -l


# ------------------------------------------------------------------------------
# 2. vcf2fas: convert per-sample GVCFs into per-gene FASTA alignments
# ------------------------------------------------------------------------------
# vcf_list.txt lists the per-sample *.concat.vcf files. The reference was split
# into chunks of 50,000 contigs to run in parallel; shown here for the whole
# reference.

ls *.concat.vcf > vcf_list.txt
vcf2fas -reference ${REFERENCE} -vcfs vcf_list.txt -gf GT -strictIUPAC 0


# ------------------------------------------------------------------------------
# 3. Filter per-gene alignments
# ------------------------------------------------------------------------------
# Run from the directory containing the .fas files. Thresholds (script defaults):
#   -g 50   max 50% missing data per site
#   -G 75   max 75% missing data per sequence
#   -H 15   max 15% heterozygous sites (IUPAC codes) per sequence
#   -l 300  min alignment length 300 bp
#   -s 30   min 30 sequences per alignment

python3 scripts/stringent_filter_fasta_alns.py -o filtered_fastas -g 50 -G 75 -H 15 -l 300 -s 30


# ------------------------------------------------------------------------------
# 4. Concatenate filtered alignments into a supermatrix
# ------------------------------------------------------------------------------
# Run from the filtered_fastas directory. Missing genes are padded with Ns.

python3 scripts/cat_alns.py -o supermatrix.fasta -t ${TAXA}


# ------------------------------------------------------------------------------
# 5. RAxML-NG: maximum likelihood tree search
# ------------------------------------------------------------------------------
# GTR+G (expanded by RAxML-NG to GTR+FO+G4m)

mpirun -np 10 raxml-ng-mpi --search \
    --msa supermatrix.fasta --msa-format FASTA \
    --model GTR+G \
    --prefix ML_search \
    --threads auto \
    --seed 12345


# ------------------------------------------------------------------------------
# 6. RAxML-NG: bootstrap replicates
# ------------------------------------------------------------------------------
# Run as ten parallel jobs (seeds 54322-54331) of 10 replicates each (100 in
# total), then combined

raxml-ng --bootstrap \
    --msa supermatrix.fasta --msa-format FASTA \
    --model GTR+G \
    --prefix bootstrap_1 \
    --threads 100 \
    --seed 54322 \
    --bs-trees 10

cat bootstrap_*.raxml.bootstraps > bootstrap_combined.raxml.bootstraps


# ------------------------------------------------------------------------------
# 7. RAxML-NG: map bootstrap support onto the ML tree
# ------------------------------------------------------------------------------
# Felsenstein bootstrap proportions (FBP) and transfer bootstrap expectation (TBE)

mpirun -np 1 raxml-ng-mpi --support \
    --tree ML_search.raxml.bestTree \
    --bs-trees bootstrap_combined.raxml.bootstraps \
    --prefix ML_support \
    --threads 4 \
    --bs-metric fbp,tbe
