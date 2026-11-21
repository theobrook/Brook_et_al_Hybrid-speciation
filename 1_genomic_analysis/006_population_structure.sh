#!/bin/bash
# ==============================================================================
# 006_population_structure.sh
# Brook et al. (2026) Not a dead end: hybrid speciation and syngameon-driven
# dynamics in an island plant radiation
#
# Model-free (PCA, MDS) and model-based (STRUCTURE, ADMIXTURE) analyses of
# population structure (Methods S2, "Model-free analysis" and "Model-based ancestry
# estimation"). Outgroups (C. prisca: C4, C43; C. repens: C77) were excluded.
#
# This file documents the commands used for each program. It is not a runnable
# pipeline.
#
# Software:
#   BCFtools v.1.21
#   PLINK v.1.90b6.21
#   ADMIXTURE v.1.3.0
#   STRUCTURE v.2.3.4
#   STRUCTURE HARVESTER
#   CLUMPP v.1.1.2
# ==============================================================================


# ------------------------------------------------------------------------------
# Variables
# ------------------------------------------------------------------------------

VCF_IN=concat_sites_stringent.vcf.gz     # linkage-filtered dataset (005_vcf_filtering.sh)


# ------------------------------------------------------------------------------
# 1. Remove outgroups
# ------------------------------------------------------------------------------

bcftools view --samples ^C43,C4,C77 -Oz -o no_outgroups.vcf.gz ${VCF_IN}


# ------------------------------------------------------------------------------
# 2. Retain biallelic SNPs and convert to PLINK binary format
# ------------------------------------------------------------------------------

bcftools view --min-alleles 2 --max-alleles 2 --types snps \
    -Oz -o biallelic_snps.vcf.gz no_outgroups.vcf.gz
bcftools index biallelic_snps.vcf.gz

plink --vcf biallelic_snps.vcf.gz --make-bed --allow-extra-chr --out snps


# ------------------------------------------------------------------------------
# 3. PCA (PLINK)
# ------------------------------------------------------------------------------

plink --vcf no_outgroups.vcf.gz --make-bed --allow-extra-chr --out snps_pca
plink --bfile snps_pca --pca 6 --allow-extra-chr --out plink_pca


# ------------------------------------------------------------------------------
# 4. MDS (PLINK)
# ------------------------------------------------------------------------------

plink --bfile snps --genome --allow-extra-chr --out plink_mds
plink --bfile snps --read-genome plink_mds.genome --allow-extra-chr \
    --cluster --mds-plot 6 --out plink_mds


# ------------------------------------------------------------------------------
# 5. STRUCTURE: prepare input
# ------------------------------------------------------------------------------
# Each diploid genotype is coded as two alleles (1/2); missing data as -9

plink --bfile snps --recode A --allow-extra-chr --out snps_raw

python3 - <<'EOF'
with open("snps_raw.raw") as f, open("structure_input.txt", "w") as out:
    next(f)
    for line in f:
        fields = line.split()
        codes = {"0": ["1", "1"], "1": ["1", "2"], "2": ["2", "2"], "NA": ["-9", "-9"]}
        genotypes = [allele for g in fields[6:] for allele in codes[g]]
        out.write(fields[1] + " " + " ".join(genotypes) + "\n")
EOF


# ------------------------------------------------------------------------------
# 6. STRUCTURE: parameter files
# ------------------------------------------------------------------------------
# Admixture model with correlated allele frequencies;
# 200,000 burn-in and 300,000 MCMC iterations

N_IND=$(wc -l < structure_input.txt)
N_LOCI=$(wc -l < snps.bim)

cat > mainparams <<EOF
#define INFILE structure_input.txt
#define OUTFILE output
#define NUMINDS ${N_IND}
#define NUMLOCI ${N_LOCI}
#define PLOIDY 2
#define MISSING -9
#define ONEROWPERIND 1
#define LABEL 1
#define POPDATA 0
#define POPFLAG 0
#define LOCDATA 0
#define PHENOTYPE 0
#define EXTRACOLS 0
#define MARKERNAMES 0
#define RECESSIVEALLELES 0
#define MAPDISTANCES 0
#define PHASED 0
#define PHASEINFO 0
#define MAXPOPS 8
#define BURNIN 200000
#define NUMREPS 300000
EOF

cat > extraparams <<EOF
#define NOADMIX 0
#define LINKAGE 0
#define USEPOPINFO 0
#define LOCPRIOR 0
#define FREQSCORR 1
#define ONEFST 0
#define INFERALPHA 1
#define ALPHA 1.0
#define POPALPHAS 0
#define UNIFPRIORALPHA 1
#define INFERLAMBDA 0
#define LAMBDA 1.0
#define COMPUTEPROB 1
#define PRINTQHAT 1
#define RANDOMIZE 0
EOF


# ------------------------------------------------------------------------------
# 7. STRUCTURE: run K = 1-8, 50 independent replicates per K
# ------------------------------------------------------------------------------
# Example for one K and one replicate; each run used a different random seed

K=2
REP=1
SEED=$((K * 100000 + REP * 1000))

structure -K ${K} -m mainparams -e extraparams \
    -i structure_input.txt \
    -o results/K${K}_rep${REP} \
    -D ${SEED}


# ------------------------------------------------------------------------------
# 8. STRUCTURE HARVESTER: Evanno delta K and CLUMPP input files
# ------------------------------------------------------------------------------

python structureHarvester.py --dir results --out harvester_results --evanno --clumpp


# ------------------------------------------------------------------------------
# 9. CLUMPP: align replicate Q-matrices (K = 2-8)
# ------------------------------------------------------------------------------
# LargeKGreedy algorithm (M = 3) with 2,000 random input orders; pairwise
# similarity statistic G' (S = 2). CLUMPP reports H', the mean pairwise G'
# across replicates (Table S1)

K=2

cat > paramfile <<EOF
DATATYPE 0
INDFILE K${K}.indfile
OUTFILE K${K}_clumpp.out
MISCFILE K${K}_clumpp.misc
K ${K}
C ${N_IND}
R 50
M 3
W 0
S 2
GREEDY_OPTION 2
REPEATS 2000
PRINT_PERMUTED_DATA 1
PERMUTED_DATAFILE K${K}_clumpp.perm_datafile
PRINT_EVERY_PERM 0
PRINT_RANDOM_INPUTORDER 0
OVERRIDE_WARNINGS 0
ORDER_BY_RUN 1
EOF

CLUMPP paramfile


# ------------------------------------------------------------------------------
# 10. ADMIXTURE: K = 2-8, 50 independent replicates per K, with cross-validation
# ------------------------------------------------------------------------------
# ADMIXTURE requires integer chromosome codes, so the contig column of the .bim
# file is set to 0. Example for one K and one replicate; each run used a
# different random seed. The K with the lowest CV error was taken as optimal.

awk '{$1="0"; print $0}' snps.bim > snps.bim.tmp && mv snps.bim.tmp snps.bim

K=2
REP=1
SEED=$((K * 100000 + REP * 1000))

admixture --cv -s ${SEED} snps.bed ${K} > K${K}_rep${REP}.log
mv snps.${K}.Q K${K}_rep${REP}.Q
mv snps.${K}.P K${K}_rep${REP}.P

# Cross-validation error
grep "CV" K${K}_rep${REP}.log


# ------------------------------------------------------------------------------
# 11. CLUMPP: align ADMIXTURE replicate Q-matrices (K = 2-8)
# ------------------------------------------------------------------------------
# Replicate .Q files are combined into a CLUMPP indfile; CLUMPP was then run
# with the same paramfile settings as for STRUCTURE (step 9)

python3 scripts/admixture_q_to_indfile.py K${K}_rep*.Q > K${K}.indfile

CLUMPP paramfile
