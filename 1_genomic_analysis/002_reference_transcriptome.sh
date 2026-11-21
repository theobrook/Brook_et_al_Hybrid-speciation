#!/bin/bash
# ==============================================================================
# 002_reference_transcriptome.sh
# Brook et al. (2026) Not a dead end: hybrid speciation and syngameon-driven
# dynamics in an island plant radiation
#
# De novo assembly of the reference transcriptome from one sample of
# C. lanceolaris (C11; Table S17) (Methods S1, "Transcriptome assembly").
#
# This file documents the commands used for each program. It is not a runnable
# pipeline.
#
# Software:
#   Trinity v.2.8.5 (including Trinity_gene_splice_modeler.py for SuperTranscripts)
#   TransDecoder v.5.7.1
#   BUSCO v.5.6.1
# ==============================================================================


# ------------------------------------------------------------------------------
# Variables
# ------------------------------------------------------------------------------

READS_DIR=results/002_trimmomatic    # quality-controlled reads for C11
ASSEMBLY_DIR=results/005_trinity
REFERENCE=Trinity.SuperTrans.fasta   # final reference used for read mapping


# ------------------------------------------------------------------------------
# 1. Trinity: de novo assembly (standard options, k-mer size 25)
# ------------------------------------------------------------------------------

Trinity --seqType fq \
    --left  ${READS_DIR}/C11_*forward_paired.fq.gz \
    --right ${READS_DIR}/C11_*reverse_paired.fq.gz \
    --max_memory 500G --CPU 100 \
    --output ${ASSEMBLY_DIR}


# ------------------------------------------------------------------------------
# 2. SuperTranscripts: collapse isoforms into one sequence per gene
# ------------------------------------------------------------------------------

Trinity_gene_splice_modeler.py --trinity_fasta ${ASSEMBLY_DIR}/Trinity.fasta

mv trinity_genes.fasta Trinity.SuperTrans.fasta
mv trinity_genes.gtf   Trinity.SuperTrans.gtf


# ------------------------------------------------------------------------------
# 3. TransDecoder: predict open reading frames (minimum 100 amino acids)
# ------------------------------------------------------------------------------
# The predicted peptides (Trinity.SuperTrans.fasta.transdecoder.pep) were used
# for functional annotation of outlier loci

TransDecoder.LongOrfs -t ${REFERENCE}
TransDecoder.Predict  -t ${REFERENCE}


# ------------------------------------------------------------------------------
# 4. BUSCO: assess completeness against three lineage datasets (Table S20)
# ------------------------------------------------------------------------------

busco --in ${REFERENCE} -m transcriptome -c 12 -l eukaryota_odb10    --out busco_eukaryota_odb10
busco --in ${REFERENCE} -m transcriptome -c 12 -l viridiplantae_odb10 --out busco_viridiplantae_odb10
busco --in ${REFERENCE} -m transcriptome -c 12 -l embryophyta_odb10  --out busco_embryophyta_odb10
