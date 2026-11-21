#!/bin/bash
# ==============================================================================
# 003_amf_identification.sh
# Brook et al. (2026) Not a dead end: hybrid speciation and syngameon-driven
# dynamics in an island plant radiation
#
# Identification of arbuscular mycorrhizal (AM) fungal ASVs (Methods S4,
# "Analysis of AM fungal community composition"):
#   1-2. BLASTn of all ASVs against the MaarjAM SSU type sequence database
#        (>= 85% identity, >= 85% query coverage, up to five hits per ASV)
#   3-5. BLASTn of every ASV with a MaarjAM hit against SILVA SSU Ref NR99, to
#        confirm that it is Glomeromycota (MaarjAM contains only AM fungi, so
#        any 18S sequence can return an AM fungal best hit)
# Hits are parsed in 004_amf_tables.R.
#
# Software:
#   BLAST+ v.2.14.1
#
# Databases:
#   MaarjAM SSU 18S rDNA type sequences, 2021 release (maarjam_database_SSU_TYPE.fasta)
#   SILVA SSU Ref NR99, release 144 (SILVA_144_SSURef_NR99_tax_silva_trunc.fasta)
# ==============================================================================


# ------------------------------------------------------------------------------
# Variables
# ------------------------------------------------------------------------------

ASVS=results/002_dada2/asvs.fasta
MAARJAM=databases/maarjam_database_SSU_TYPE.fasta
SILVA=databases/SILVA_144_SSURef_NR99_tax_silva_trunc.fasta
OUT_DIR=results/003_blast


# ------------------------------------------------------------------------------
# 1. MaarjAM BLAST database
# ------------------------------------------------------------------------------

makeblastdb \
    -in ${MAARJAM} \
    -dbtype nucl \
    -out databases/db_maarjam \
    -title "MaarjAM_SSU_TYPE"


# ------------------------------------------------------------------------------
# 2. BLASTn against MaarjAM
# ------------------------------------------------------------------------------

blastn \
    -query ${ASVS} \
    -db databases/db_maarjam \
    -out ${OUT_DIR}/blast_maarjam.txt \
    -outfmt "6 qseqid sseqid pident length qlen qcovs evalue bitscore stitle" \
    -perc_identity 85 \
    -qcov_hsp_perc 85 \
    -max_target_seqs 5 \
    -num_threads 8


# ------------------------------------------------------------------------------
# 3. SILVA BLAST database
# ------------------------------------------------------------------------------
# SILVA sequences are RNA: convert U to T on sequence lines first

sed '/^>/!y/Uu/Tt/' ${SILVA} > databases/silva_ssu_nr99_dna.fasta

makeblastdb \
    -in databases/silva_ssu_nr99_dna.fasta \
    -dbtype nucl \
    -out databases/db_silva \
    -title "SILVA_SSURef_NR99"


# ------------------------------------------------------------------------------
# 4. Query = ASVs with a MaarjAM hit
# ------------------------------------------------------------------------------

cut -f1 ${OUT_DIR}/blast_maarjam.txt | sort -u > ${OUT_DIR}/maarjam_hit_asvs.txt

awk 'NR == FNR {keep[">"$1]; next} /^>/ {p = ($1 in keep)} p' \
    ${OUT_DIR}/maarjam_hit_asvs.txt ${ASVS} > ${OUT_DIR}/silva_query.fasta


# ------------------------------------------------------------------------------
# 5. BLASTn against SILVA
# ------------------------------------------------------------------------------

blastn \
    -query ${OUT_DIR}/silva_query.fasta \
    -db databases/db_silva \
    -out ${OUT_DIR}/blast_silva.txt \
    -outfmt "6 qseqid sseqid pident length qlen qcovs evalue bitscore stitle" \
    -max_target_seqs 20 \
    -evalue 1e-10 \
    -num_threads 8
