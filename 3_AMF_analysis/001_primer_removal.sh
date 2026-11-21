#!/bin/bash
# ==============================================================================
# 001_primer_removal.sh
# Brook et al. (2026) Not a dead end: hybrid speciation and syngameon-driven
# dynamics in an island plant radiation
#
# Removal of the 18S V4 primers (and the heterogeneity spacers 5' of them)
# from paired root and soil metabarcoding reads (Methods S4, "Analysis of AM
# fungal community composition").
#
# This file documents the commands used. It is not a runnable pipeline: the
# command is shown for one example sample and was run for every sample.
#
# Software:
#   cutadapt v.5.2
# ==============================================================================


# ------------------------------------------------------------------------------
# Variables
# ------------------------------------------------------------------------------

RAW_DIR=data/18s_raw                 # Biomarker Technologies FASTQ files
OUT_DIR=results/001_cutadapt

FWD=CCAGCASCYGCGGTAATTCC             # TAReuk454FWD1
REV=ACTTTCGTTCTTGATYRA               # TAReukREV3

SAMPLE=Unknown_CA007_001M0001        # example sample


# ------------------------------------------------------------------------------
# 1. cutadapt
# ------------------------------------------------------------------------------
# Non-anchored -g/-G remove each primer and everything 5' of it (spacers).
# Maximum error rate 0.1, minimum primer overlap 15 bp. Pairs in which either
# primer was not found, or either read is < 100 bp after trimming, are discarded.

cutadapt \
    -g ${FWD} \
    -G ${REV} \
    -e 0.1 \
    --overlap 15 \
    --discard-untrimmed \
    --pair-filter any \
    --minimum-length 100 \
    -o ${OUT_DIR}/${SAMPLE}_1.fq.gz \
    -p ${OUT_DIR}/${SAMPLE}_2.fq.gz \
    ${RAW_DIR}/${SAMPLE}_1.fq.gz \
    ${RAW_DIR}/${SAMPLE}_2.fq.gz \
    > ${OUT_DIR}/${SAMPLE}.log
