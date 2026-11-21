#!/bin/bash
# ==============================================================================
# 003_maxent.sh
# Brook et al. (2026) Not a dead end: hybrid speciation and syngameon-driven
# dynamics in an island plant radiation
#
# Species distribution models (Methods S3, "Species distribution analysis";
# Tables S6, S9 and S10; Fig. S5).
#
# The models in the paper were run in the MaxEnt graphical interface. These
# commands reproduce the same settings from the command line:
#   - full (unthinned) occurrence dataset; MaxEnt's default removal of
#     duplicate presence records within a grid cell
#   - nine environmental layers (10 x 10 m), geology and vegetation categorical
#   - 10-fold cross-validation (5-fold for C. inopinata, run separately)
#   - convergence threshold 10^-5, maximum 500 iterations
#   - cloglog output, response curves and jackknife tests of variable importance
#   - no bias file
# randomseed=true was used, so replicate splits (and therefore exact values)
# differ between runs.
#
# Software:
#   MaxEnt v.3.4 (maxent.jar), Java
# ==============================================================================


# ------------------------------------------------------------------------------
# Variables
# ------------------------------------------------------------------------------

MAXENT=maxent/maxent.jar
SAMPLES=data/species_occurrence.csv                              # Species, x, y (full dataset)
LAYERS=data/ascii_layers                                         # nine .asc layers
OUT_DIR=maxent/no_thinning
OUT_DIR_INOP=maxent/no_thinning_inopinata

mkdir -p ${OUT_DIR} ${OUT_DIR_INOP}


# ------------------------------------------------------------------------------
# 1. C. huttoniana, C. lanceolaris, C. ptotopetra, C. putida, C. savolainenii
#    (10-fold cross-validation)
# ------------------------------------------------------------------------------
# togglespeciesselected deselects species whose names begin with the prefix:
# C. inopinata is run separately below and C. prisca was not modelled

java -mx4g -jar ${MAXENT} \
    environmentallayers=${LAYERS} \
    samplesfile=${SAMPLES} \
    outputdirectory=${OUT_DIR} \
    togglelayertype=geology \
    togglelayertype=vegetation \
    togglespeciesselected=inopinata \
    togglespeciesselected=prisca \
    replicates=10 \
    replicatetype=crossvalidate \
    randomseed=true \
    maximumiterations=500 \
    convergencethreshold=0.00001 \
    outputformat=cloglog \
    responsecurves=true \
    jackknife=true \
    writeplotdata=true \
    redoifexists \
    autorun \
    visible=false


# ------------------------------------------------------------------------------
# 2. C. inopinata (5-fold cross-validation, small sample size)
# ------------------------------------------------------------------------------
# Samples file with C. inopinata records only

awk -F',' 'NR == 1 || $1 ~ /inopinata/' ${SAMPLES} > ${OUT_DIR_INOP}/inopinata_samples.csv

java -mx4g -jar ${MAXENT} \
    environmentallayers=${LAYERS} \
    samplesfile=${OUT_DIR_INOP}/inopinata_samples.csv \
    outputdirectory=${OUT_DIR_INOP} \
    togglelayertype=geology \
    togglelayertype=vegetation \
    replicates=5 \
    replicatetype=crossvalidate \
    randomseed=true \
    maximumiterations=500 \
    convergencethreshold=0.00001 \
    outputformat=cloglog \
    responsecurves=true \
    jackknife=true \
    writeplotdata=true \
    redoifexists \
    autorun \
    visible=false


# ------------------------------------------------------------------------------
# Outputs used downstream
# ------------------------------------------------------------------------------
# <species>_avg.asc                 habitat suitability averaged across
#                                   replicates (Fig. S5; Schoener's D and
#                                   Warren's I, Table S9)
# maxentResults.csv                 test AUC per replicate and variable
#                                   contributions (Tables S6 and S10)
