#!/bin/bash
# ==============================================================================
# 003_read_mapping.sh
# Brook et al. (2026) Not a dead end: hybrid speciation and syngameon-driven
# dynamics in an island plant radiation
#
# Mapping reads to the reference transcriptome and preparing BAM files for
# variant calling (Methods S1, "Mapping, alignment, and variant calling").
#
# This file documents the commands used for each program. It is not a runnable
# pipeline: per-sample commands are shown for one example sample and were run
# for every sample.
#
# Software:
#   STAR v.2.7.11b
#   Picard v.2.18.21
#   GATK v.4.5.0.0
#   SAMtools v.1.21
# ==============================================================================


# ------------------------------------------------------------------------------
# Variables
# ------------------------------------------------------------------------------

REFERENCE=Trinity.SuperTrans.fasta
DICT=Trinity.SuperTrans.dict
READS_DIR=results/003_concatenated
GENOME_DIR=results/006_star_index
STAR_DIR=results/007_star
RG_DIR=results/008_picard_readgroups
DEDUP_DIR=results/009_picard_dedup
SPLIT_DIR=results/010_gatk_splitncigar

SAMPLE=C10                           # example sample


# ------------------------------------------------------------------------------
# 1. STAR: build the reference index
# ------------------------------------------------------------------------------

STAR --runMode genomeGenerate \
    --runThreadN 90 \
    --limitGenomeGenerateRAM 250000000000 \
    --genomeDir ${GENOME_DIR} \
    --genomeFastaFiles ${REFERENCE} \
    --genomeSAindexNbases 11


# ------------------------------------------------------------------------------
# 2. STAR: map each sample separately
# ------------------------------------------------------------------------------

STAR --runThreadN 50 \
    --genomeDir ${GENOME_DIR} \
    --readFilesIn ${READS_DIR}/${SAMPLE}_forward_paired.fq.gz \
                  ${READS_DIR}/${SAMPLE}_reverse_paired.fq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --outFileNamePrefix ${STAR_DIR}/${SAMPLE}_


# ------------------------------------------------------------------------------
# 3. Picard: add read groups
# ------------------------------------------------------------------------------

picard -Xmx5g AddOrReplaceReadGroups \
    I=${STAR_DIR}/${SAMPLE}_Aligned.sortedByCoord.out.bam \
    O=${RG_DIR}/${SAMPLE}_rg.bam \
    RGID=4 RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=${SAMPLE}


# ------------------------------------------------------------------------------
# 4. Picard: remove PCR and optical duplicates
# ------------------------------------------------------------------------------

picard -Xmx5g MarkDuplicates \
    I=${RG_DIR}/${SAMPLE}_rg.bam \
    O=${DEDUP_DIR}/${SAMPLE}_processed.bam \
    M=${DEDUP_DIR}/${SAMPLE}_metrics.txt \
    REMOVE_DUPLICATES=true


# ------------------------------------------------------------------------------
# 5. Reference dictionary (GATK) and FASTA index (SAMtools)
# ------------------------------------------------------------------------------

gatk CreateSequenceDictionary -R ${REFERENCE}
samtools faidx ${REFERENCE}


# ------------------------------------------------------------------------------
# 6. GATK SplitNCigarReads: split reads spanning splice junctions and
#    hard-clip overhangs
# ------------------------------------------------------------------------------

gatk --java-options "-Xmx80g" SplitNCigarReads \
    -R ${REFERENCE} \
    --sequence-dictionary ${DICT} \
    -I ${DEDUP_DIR}/${SAMPLE}_processed.bam \
    -O ${SPLIT_DIR}/${SAMPLE}_s.bam
