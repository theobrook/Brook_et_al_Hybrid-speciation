#!/bin/bash
# ==============================================================================
# 004_variant_calling.sh
# Brook et al. (2026) Not a dead end: hybrid speciation and syngameon-driven
# dynamics in an island plant radiation
#
# Per-sample variant calling, merging and removal of high-coverage sites
# (Methods S1, "Mapping, alignment, and variant calling").
#
# This file documents the commands used for each program. It is not a runnable
# pipeline: per-sample commands are shown for one example sample and were run
# for every sample.
#
# Software:
#   BCFtools v.1.21
# ==============================================================================


# ------------------------------------------------------------------------------
# Variables
# ------------------------------------------------------------------------------

REFERENCE=Trinity.SuperTrans.fasta
BAM_DIR=results/010_gatk_splitncigar
BCF_DIR=results/011_bcftools_call

SAMPLE=C10                           # example sample


# ------------------------------------------------------------------------------
# 1. Call variants per sample (GVCF output, max per-sample depth 300)
# ------------------------------------------------------------------------------
# mpileup -> call (multiallelic caller, GQ annotation) -> norm (split
# multiallelic sites, left-align indels) -> fill-tags -> compressed BCF

bcftools mpileup -Ou --gvcf 20000 -d 300 --fasta-ref ${REFERENCE} ${BAM_DIR}/${SAMPLE}_s.bam | \
bcftools call -Ou -m -f GQ --gvcf 20000 | \
bcftools norm -m+any --fasta-ref ${REFERENCE} | \
bcftools +fill-tags -- -t all | \
bcftools view -Ob -o ${BCF_DIR}/${SAMPLE}call.bcf

bcftools index ${BCF_DIR}/${SAMPLE}call.bcf


# ------------------------------------------------------------------------------
# 2. Merge all samples into one file
# ------------------------------------------------------------------------------

ls ${BCF_DIR}/*.bcf > bcf_file_list.txt
bcftools merge --file-list bcf_file_list.txt -Ob -o merged.bcf


# ------------------------------------------------------------------------------
# 3. Add INFO tags to the merged file and remove sites with excessive coverage
#    (INFO/DP >= 1260; potential collapsed paralogues)
# ------------------------------------------------------------------------------

bcftools +fill-tags merged.bcf -Ob -o merged_tags.bcf -- -t all
bcftools view -i 'INFO/DP < 1260' merged_tags.bcf -Ob -o merged_tags.MaxDPfilt.bcf
