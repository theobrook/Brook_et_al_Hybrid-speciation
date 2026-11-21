# Brook et al. (2026) Not a dead end: hybrid speciation and syngameon-driven dynamics in an island plant radiation

This repository contains the analysis scripts for Brook et al. (2026), a study of hybrid speciation and syngameon dynamics in the radiation of *Coprosma* on Lord Howe Island. The scripts cover the genomic analyses (RNA-seq), the aboveground ecological analyses and the analyses of root arbuscular mycorrhizal (AM) fungal communities (18S metabarcoding). Methods S1–S4 in the Supplementary Information describe each step in full.

The bash (`.sh`) files record the commands and settings used for each program. They are not runnable pipelines: each command is shown for one example sample, K or species pair, and the analyses were run on an HPC cluster. The R scripts can be run once their input paths are set in the `Settings` section at the top of each file.

## 1_genomic_analysis

RNA-seq read processing, reference transcriptome, variant calling and population genomic analyses (Methods S1 and S2).

1. `001_read_qc.sh`: adapter and quality trimming (Trim Galore, Trimmomatic) and read quality checks (FastQC).
2. `002_reference_transcriptome.sh`: de novo assembly of the reference transcriptome (Trinity), SuperTranscripts, ORF prediction (TransDecoder) and completeness (BUSCO).
3. `003_read_mapping.sh`: read mapping to the reference (STAR), duplicate marking (Picard) and splitting of reads spanning introns (GATK SplitNCigarReads).
4. `004_variant_calling.sh`: per-sample variant calling (BCFtools mpileup/call), merging and removal of high-coverage sites.
5. `005_vcf_filtering.sh`: variant filtering (VCFtools) and linkage filtering (one SNP per contig).
6. `006_population_structure.sh`: PCA and MDS (PLINK), STRUCTURE, STRUCTURE HARVESTER, CLUMPP and ADMIXTURE.
7. `006_structure_plots.R`: ADMIXTURE and STRUCTURE bar plots of the best replicate per K (Fig. 1D–E).
8. `006_triangle_plots.R`: triangle plots of hybrid index against interclass heterozygosity (triangulaR; Fig. 2).
9. `007_genetic_diversity.sh`: observed and expected heterozygosity and F<sub>IS</sub> (PLINK), all-sites VCF, π and D<sub>XY</sub> (pixy), and pairwise weighted F<sub>ST</sub> (VCFtools).
10. `007_pixy_summary.R`: genome-wide π per species and D<sub>XY</sub> per species pair (Tables S2 and S3).
11. `007_private_alleles.R`: private alleles and fixed private differences per species, with rarefaction.
12. `007_fst_hybrid_test.R`: phylogenetically corrected F<sub>ST</sub> and the Wilcoxon test of hybrid–parent versus non-hybrid pairs.
13. `008_phylogenetics.sh`: maximum likelihood phylogeny and bootstrap support (RAxML-NG).
14. `008_ploidy_multiallelic.sh`: counts of all genotyped and multiallelic sites per individual in the merged BCF (BCFtools; Table S16).
15. `008_ploidy_multiallelic.R`: proportion of multiallelic sites per individual and Welch's t-test between diploid species and species of unknown ploidy (Methods S6; Results S3).
16. `009_selection_outliers.sh`: functional annotation (eggNOG-mapper), per-site F<sub>ST</sub> (pixy), BLASTp against Araport11 and keyword screening of outliers.
17. `009_fst_outliers.R`: top 10% and top 5% F<sub>ST</sub> outliers per species pair (Table S12).
18. `data/multiallelic_sites.csv`: callable and multiallelic sites per individual (Table S16), used in `008_ploidy_multiallelic.R`.
19. `scripts/`: Python helper scripts called by the files above:
    - `admixture_q_to_indfile.py`: converts ADMIXTURE Q files to CLUMPP input.
    - `cat_alns.py` and `stringent_filter_fasta_alns.py`: build and filter the alignment for RAxML-NG.
    - `annotate_outliers_eggnog.py`, `tair_annotation.py` and `merge_blast_annotations.py`: annotate outlier transcripts with eggNOG-mapper and TAIR.
    - `summarise_outliers.py`: summarises outliers and flags local adaptation and reproductive isolation candidates (Tables S13 and S14).

## 2_ecological_analysis

Aboveground ecological differentiation and its relationship with reproductive isolation (Methods S3).

1. `001_spatial_thinning.R`: spatially thinned occurrence dataset (one record per species per 10 x 10 m cell; Table S18).
2. `scripts/extract_values.py`: extracts the values of the nine environmental layers at each occurrence point.
3. `002_environmental_pca.R`: PCA of the nine environmental variables at the spatially thinned occurrence points (Table S5; Figs 3 and S4), and Euclidean distances between species centroids across all nine PCs.
4. `002_environmental_differentiation.R`: global and pairwise PERMANOVAs of environmental differentiation (Table S7) and pairwise Wilcoxon tests on each PC (Table S8).
5. `003_maxent.sh`: species distribution models (MaxEnt), written as command-line equivalents of the settings used in the MaxEnt graphical interface.
6. `003_niche_overlap.R`: Schoener's *D* and Warren's *I* between all species pairs from the averaged MaxEnt suitability surfaces (Table S9).
7. `004_fst_niche_correlations.R`: Pearson correlations between phylogenetically corrected F<sub>ST</sub> and centroid distance, Schoener's *D* and Warren's *I*.
8. `005_flowering_time.R`: Kruskal–Wallis and pairwise Wilcoxon tests of flowering time (Table S15), and flowering-time plots (Figs S8 and S9).
9. `data/species_occurrence.csv`: occurrence records for each species (full dataset; Table S18), in projected coordinates (metres), used in `001_spatial_thinning.R` and `003_maxent.sh`.
10. `data/ascii_layers/`: the nine environmental layers (10 x 10 m, ASCII grid format) from Papadopulos et al. (2013), used in `scripts/extract_values.py` and `003_maxent.sh`.
11. `data/flowering_time.csv`: flowering-time records from herbarium specimens (NSW and K), field collections and field observations, used in `005_flowering_time.R`.

The species distribution models were run in the MaxEnt v.3.4 graphical interface; `003_maxent.sh` gives the equivalent command-line settings.

## 3_AMF_analysis

Root AM fungal community composition from 18S V4 metabarcoding (Methods S4).

1. `001_primer_removal.sh`: removal of primers and heterogeneity spacers (cutadapt).
2. `002_dada2_asvs.R`: quality filtering, denoising, read merging and chimera removal (DADA2), and rarefaction of all ASVs (Table S19).
3. `003_amf_identification.sh`: BLASTn of ASVs against MaarjAM, and of MaarjAM hits against SILVA.
4. `004_amf_tables.R`: AM fungal ASV and virtual taxon (VT) tables from the SILVA-confirmed MaarjAM hits.
5. `005_amf_composition.R`: global PERMANOVAs (Table S11), PCA (Fig. S7) and genus composition (Fig. S6).
6. `006_qpcr_colonisation.R`: assumption checks, Kruskal–Wallis and Dunn's tests of root AM fungal abundance from qPCR (Results S1).
7. `data/qpcr_results.csv`: qPCR Ct values and AM fungal DNA estimates per root and soil sample (Table S22), used in `006_qpcr_colonisation.R`.

## Data

- Code and ecological data: https://github.com/theobrook/Brook_et_al_Hybrid-speciation
- RNA-seq reads: NCBI BioProject [PRJNA1560389](https://www.ncbi.nlm.nih.gov/bioproject/1560389)
- 18S metabarcoding reads: NCBI BioProject [PRJNA1560391](https://www.ncbi.nlm.nih.gov/bioproject/1560391)
- Environmental layers: Papadopulos et al. (2013) *Journal of Evolutionary Biology* 26, 733–745.
- Occurrence records and herbarium vouchers: Tables S17, S18 and S21.
