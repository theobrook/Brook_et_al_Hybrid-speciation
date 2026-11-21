#!/bin/bash
# ==============================================================================
# 009_selection_outliers.sh
# Brook et al. (2026) Not a dead end: hybrid speciation and syngameon-driven
# dynamics in an island plant radiation
#
# Per-site FST outliers between species pairs and their functional annotation
# (main Methods, "Selection analyses"). Individuals with substantial mixed
# ancestry at K = 6 (C67, C70) were excluded; outgroups were excluded when
# identifying outliers (15 pairwise comparisons).
#
# This file documents the commands used for each program. It is not a runnable
# pipeline: pixy is shown for one example species pair and was run for every
# pair.
#
# Software:
#   eggNOG-mapper v.2.1.13 (eggNOG database v.5.0.2)
#   pixy v.2.0.0.beta14
#   R (009_fst_outliers.R)
#   BLAST+ v.2.15.0
#   Python 3 (scripts/)
# ==============================================================================


# ------------------------------------------------------------------------------
# Variables
# ------------------------------------------------------------------------------

VCF=allsites_stringent.vcf.gz          # biallelic SNPs + invariant sites (007_genetic_diversity.sh)
POPS_FILE=populations.txt              # two columns (sample, species), excluding C67 and C70
PEPTIDES=Trinity.SuperTrans.fasta.transdecoder.pep   # predicted peptides (002_reference_transcriptome.sh)
ARAB_PEP=Araport11_pep_20250411        # Araport11 protein sequences, all isoforms, release 2025-04-11 (TAIR)
ARAB_DB=arabidopsis_araport11_2025     # BLAST database built from ${ARAB_PEP}
OUT=results/outliers

POP1=huttoniana                        # example species pair
POP2=ptotopetra


# ------------------------------------------------------------------------------
# 0. eggNOG-mapper: functional annotation of the reference transcriptome
# ------------------------------------------------------------------------------
# Database files (eggNOG v5.0.2) downloaded into results/eggnog from
# http://eggnog5.embl.de/download/emapperdb-5.0.2/ :
#   eggnog.db.gz, eggnog_proteins.dmnd.gz, eggnog.taxa.tar.gz

emapper.py \
    -i ${PEPTIDES} \
    --itype proteins \
    -o coprosma_annotation \
    --data_dir results/eggnog \
    --cpu 4


# ------------------------------------------------------------------------------
# 1. pixy: per-site Weir-Cockerham FST for each species pair
# ------------------------------------------------------------------------------

awk -v p1=${POP1} -v p2=${POP2} '$2 == p1 || $2 == p2' ${POPS_FILE} > ${POP1}_${POP2}_pops.txt

pixy --stats fst \
    --vcf ${VCF} \
    --populations ${POP1}_${POP2}_pops.txt \
    --window_size 1 \
    --output_folder ${OUT}/fst \
    --output_prefix ${POP1}_vs_${POP2} \
    --n_cores 4


# ------------------------------------------------------------------------------
# 2. Identify top 10% and top 5% FST outliers per pair, and map outlier sites
#    to transcripts (R)
# ------------------------------------------------------------------------------
# Ties at each threshold are retained (quantile type = 1)

Rscript 009_fst_outliers.R


# ------------------------------------------------------------------------------
# 3. Annotate outlier transcripts with eggNOG-mapper output and flag keywords
#    for local adaptation and RI (Table S13)
# ------------------------------------------------------------------------------

python3 scripts/annotate_outliers_eggnog.py


# ------------------------------------------------------------------------------
# 4. BLASTp of outlier proteins against the Araport11 proteome
# ------------------------------------------------------------------------------

# 4a. Local protein BLAST database from the Araport11 proteome

gunzip -k ${ARAB_PEP}.gz

makeblastdb \
    -in ${ARAB_PEP} \
    -dbtype prot \
    -parse_seqids \
    -out ${ARAB_DB} \
    -title "Araport11_proteome"

# 4b. Peptide sequences of all unique outlier transcripts: the union of the
#     top 10% lists of all 15 pairs (009_fst_outliers.R), which includes the
#     top 5% outliers

cat ${OUT}/003_annotation/pairwise_outliers/*/FST_top10/outlier_transcripts.txt \
    | sort -u > ${OUT}/unique_outlier_transcripts.txt

awk 'NR == FNR {keep[">"$1]; next} /^>/ {p = ($1 in keep)} p' \
    ${OUT}/unique_outlier_transcripts.txt ${PEPTIDES} > outlier_proteins.fasta

# 4c. BLASTp

blastp \
    -query outlier_proteins.fasta \
    -db ${ARAB_DB} \
    -out ${OUT}/006_blast_arabidopsis/blast_results.tsv \
    -evalue 1e-5 \
    -max_target_seqs 5 \
    -num_threads 8 \
    -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen stitle"

# Best hit per query = highest bitscore; query coverage = length / qlen x 100
# (written to 006_blast_arabidopsis/best_arabidopsis_hits.tsv)


# ------------------------------------------------------------------------------
# 5. Merge eggNOG and BLAST annotations, add Araport11 functional descriptions,
#    and summarise outliers across the 15 pairwise comparisons
# ------------------------------------------------------------------------------

python3 scripts/merge_blast_annotations.py
python3 scripts/tair_annotation.py
python3 scripts/summarise_outliers.py
