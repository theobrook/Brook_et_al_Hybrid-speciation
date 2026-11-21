#!/bin/bash
# ==============================================================================
# 001_read_qc.sh
# Brook et al. (2026) Not a dead end: hybrid speciation and syngameon-driven
# dynamics in an island plant radiation
#
# Quality control of RNA-seq reads (Methods S1, "Quality control").
# Reads were delivered by Novogene UK after their initial filtering.
#
# This file documents the commands used for each program. It is not a runnable
# pipeline: commands are shown for one example sample and were run for every
# sample (and every sequencing lane) on an HPC cluster.
#
# Software:
#   Trim Galore v.0.6.10
#   Trimmomatic v.0.39
#   FastQC v.0.12.1
# ==============================================================================


# ------------------------------------------------------------------------------
# Variables
# ------------------------------------------------------------------------------

RAW_DIR=data/rnaseq_raw              # Novogene FASTQ files (one folder per sample)
TRIMGALORE_DIR=results/001_trimgalore
TRIMMOMATIC_DIR=results/002_trimmomatic
CONCAT_DIR=results/003_concatenated
FASTQC_DIR=results/004_fastqc
ADAPTERS=Trimmomatic-0.39/adapters/TruSeq3-PE-2.fa

SAMPLE=C10                           # example sample
LANE=L3                              # example lane


# ------------------------------------------------------------------------------
# 1. Trim Galore: remove remaining Illumina adapters and quality-trim
# ------------------------------------------------------------------------------
# Defaults: adapter overlap >= 1 bp, 3' quality trimming at Phred 20,
# discard reads < 20 bp

trim_galore --illumina --paired --cores 35 \
    ${RAW_DIR}/${SAMPLE}/${SAMPLE}_*_${LANE}_1.fq.gz \
    ${RAW_DIR}/${SAMPLE}/${SAMPLE}_*_${LANE}_2.fq.gz \
    -o ${TRIMGALORE_DIR}


# ------------------------------------------------------------------------------
# 2. Trimmomatic: second round of adapter and quality trimming
# ------------------------------------------------------------------------------

trimmomatic PE -threads 8 -phred33 \
    ${TRIMGALORE_DIR}/${SAMPLE}_*_${LANE}_1_val_1.fq.gz \
    ${TRIMGALORE_DIR}/${SAMPLE}_*_${LANE}_2_val_2.fq.gz \
    ${TRIMMOMATIC_DIR}/${SAMPLE}_${LANE}_forward_paired.fq.gz \
    ${TRIMMOMATIC_DIR}/${SAMPLE}_${LANE}_forward_unpaired.fq.gz \
    ${TRIMMOMATIC_DIR}/${SAMPLE}_${LANE}_reverse_paired.fq.gz \
    ${TRIMMOMATIC_DIR}/${SAMPLE}_${LANE}_reverse_unpaired.fq.gz \
    ILLUMINACLIP:${ADAPTERS}:2:30:10 \
    LEADING:5 TRAILING:5 SLIDINGWINDOW:4:5 MINLEN:75 TOPHRED33


# ------------------------------------------------------------------------------
# 3. Combine reads from samples sequenced across more than one lane/run
# ------------------------------------------------------------------------------
# Only paired reads were used downstream

cat ${TRIMMOMATIC_DIR}/${SAMPLE}_*_forward_paired.fq.gz > ${CONCAT_DIR}/${SAMPLE}_forward_paired.fq.gz
cat ${TRIMMOMATIC_DIR}/${SAMPLE}_*_reverse_paired.fq.gz > ${CONCAT_DIR}/${SAMPLE}_reverse_paired.fq.gz


# ------------------------------------------------------------------------------
# 4. FastQC: check final read quality
# ------------------------------------------------------------------------------

fastqc -t 8 -o ${FASTQC_DIR} ${CONCAT_DIR}/*_paired.fq.gz
