################################################################################
# 002_dada2_asvs.R
# Brook et al. (2026) Not a dead end: hybrid speciation and syngameon-driven
# dynamics in an island plant radiation
#
# Primer-trimmed 18S reads (001_primer_removal.sh) -> amplicon sequence
# variants (ASVs) and rarefied ASV tables (Methods S4, "Analysis of AM fungal
# community composition"; Table S19).
#
#   1. Quality profiles of trimmed reads, per sequencing run
#   2. filterAndTrim: truncate to 210 / 200 bp; maxN = 0; maxEE = 2; PhiX removed
#   3. Error models learned separately for each NovaSeq run, with a loess fit
#      adapted to binned quality scores (dada2 GitHub issue #1307)
#   4. Sample-by-sample denoising, read merging (minimum overlap 12 bp)
#   5. Sequence tables of the two runs combined; ASVs > 30 bp from the modal
#      length removed; chimeras removed (consensus method)
#   6. ASV names, FASTA and ASV table
#   7. All ASVs rarefied (rrarefy) separately for roots and soils, to the
#      minimum depth within each sample type and dataset (batch 1; batches 1 + 2)
#
# Software: R v.4.x, dada2 v.1.30.0, vegan v.2.7.1
################################################################################

library(dada2)
library(ShortRead)
library(Biostrings)
library(tidyverse)
library(vegan)


# ------------------------------------------------------------------------------
# Settings
# ------------------------------------------------------------------------------

cut_dir      <- "results/001_cutadapt"       # *_1.fq.gz / *_2.fq.gz
sample_sheet <- "sample_sheet.csv"           # seqtab_name, Sample_ID, Species, Sample_type, Batch
out_dir      <- "results/002_dada2"
dir.create(file.path(out_dir, "filtered"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "rarefied"), recursive = TRUE, showWarnings = FALSE)

truncF        <- 210
truncR        <- 200
min_overlap   <- 12
len_tolerance <- 30       # keep ASVs within modal length +/- 30 bp
rare_seed     <- 123
n_threads     <- TRUE


# ------------------------------------------------------------------------------
# Files and sequencing runs
# ------------------------------------------------------------------------------

fnFs <- sort(list.files(cut_dir, "_1\\.fq\\.gz$", full.names = TRUE))
fnRs <- sort(list.files(cut_dir, "_2\\.fq\\.gz$", full.names = TRUE))
samples <- sub("_1\\.fq\\.gz$", "", basename(fnFs))
names(fnFs) <- names(fnRs) <- samples

# Sequencing run = instrument:run:flowcell from the first read header
get_run <- function(f) {
  con <- gzfile(f); on.exit(close(con))
  paste(strsplit(sub("^@", "", readLines(con, n = 1)), ":")[[1]][1:3], collapse = ":")
}
run_id <- vapply(fnFs, get_run, character(1))
names(run_id) <- samples


# ------------------------------------------------------------------------------
# 1. Quality profiles (inspected to choose the truncation lengths)
# ------------------------------------------------------------------------------

for (r in unique(run_id)) {
  s   <- samples[run_id == r]
  tag <- gsub("[^A-Za-z0-9]", "_", r)
  ggsave(file.path(out_dir, paste0("quality_F_", tag, ".png")),
         plotQualityProfile(fnFs[s], aggregate = TRUE), width = 10, height = 7)
  ggsave(file.path(out_dir, paste0("quality_R_", tag, ".png")),
         plotQualityProfile(fnRs[s], aggregate = TRUE), width = 10, height = 7)
}


# ------------------------------------------------------------------------------
# 2. Filter and trim
# ------------------------------------------------------------------------------

filtFs <- file.path(out_dir, "filtered", paste0(samples, "_F_filt.fastq.gz"))
filtRs <- file.path(out_dir, "filtered", paste0(samples, "_R_filt.fastq.gz"))
names(filtFs) <- names(filtRs) <- samples

filt <- filterAndTrim(fnFs, filtFs, fnRs, filtRs,
                      truncLen = c(truncF, truncR), maxN = 0, maxEE = c(2, 2),
                      truncQ = 2, rm.phix = TRUE, compress = TRUE,
                      multithread = n_threads)
rownames(filt) <- samples
kept <- samples[filt[, "reads.out"] > 0]


# ------------------------------------------------------------------------------
# 3. Error model for NovaSeq binned quality scores
# ------------------------------------------------------------------------------
# loess weighted by total counts, degree 1, and forced to be non-increasing
# with quality score

loess_errfun_binned <- function(trans) {
  qq  <- as.numeric(colnames(trans))
  est <- matrix(0, nrow = 0, ncol = length(qq))
  for (nti in c("A", "C", "G", "T")) for (ntj in c("A", "C", "G", "T")) if (nti != ntj) {
    errs  <- trans[paste0(nti, "2", ntj), ]
    tot   <- colSums(trans[paste0(nti, "2", c("A", "C", "G", "T")), ])
    rlogp <- log10((errs + 1) / tot)
    rlogp[is.infinite(rlogp)] <- NA
    df    <- data.frame(q = qq, errs = errs, tot = tot, rlogp = rlogp)
    pred  <- predict(loess(rlogp ~ q, df, weights = log10(tot), span = 2, degree = 1), qq)
    hi <- max(which(!is.na(pred))); lo <- min(which(!is.na(pred)))
    pred[seq_along(pred) > hi] <- pred[[hi]]
    pred[seq_along(pred) < lo] <- pred[[lo]]
    est <- rbind(est, 10^pred)
  }
  est[est > 0.25] <- 0.25
  est[est < 1e-7] <- 1e-7
  est <- t(apply(est, 1, function(x) rev(cummax(rev(x)))))
  err <- rbind(1 - colSums(est[1:3, ]), est[1:3, ],
               est[4, ], 1 - colSums(est[4:6, ]), est[5:6, ],
               est[7:8, ], 1 - colSums(est[7:9, ]), est[9, ],
               est[10:12, ], 1 - colSums(est[10:12, ]))
  rownames(err) <- paste0(rep(c("A", "C", "G", "T"), each = 4), "2", c("A", "C", "G", "T"))
  colnames(err) <- colnames(trans)
  err
}


# ------------------------------------------------------------------------------
# 4. Denoise and merge, separately for each sequencing run
# ------------------------------------------------------------------------------

get_n   <- function(x) sum(getUniques(x))
seqtabs <- list()
counts  <- list()

for (r in unique(run_id[kept])) {
  s <- kept[run_id[kept] == r]

  errF <- learnErrors(filtFs[s], errorEstimationFunction = loess_errfun_binned,
                      nbases = 1e8, randomize = TRUE, multithread = n_threads)
  errR <- learnErrors(filtRs[s], errorEstimationFunction = loess_errfun_binned,
                      nbases = 1e8, randomize = TRUE, multithread = n_threads)

  dF <- dada(filtFs[s], err = errF, pool = FALSE, multithread = n_threads)
  dR <- dada(filtRs[s], err = errR, pool = FALSE, multithread = n_threads)
  mg <- mergePairs(dF, filtFs[s], dR, filtRs[s], minOverlap = min_overlap)

  seqtabs[[r]] <- makeSequenceTable(mg)
  counts[[r]]  <- tibble(sample = s,
                         denoisedF = sapply(dF, get_n),
                         denoisedR = sapply(dR, get_n),
                         merged    = sapply(mg, get_n))
}


# ------------------------------------------------------------------------------
# 5. Combine runs, length filter, chimera removal
# ------------------------------------------------------------------------------

seqtab <- mergeSequenceTables(tables = unname(seqtabs))

len_mode <- as.numeric(names(which.max(table(nchar(getSequences(seqtab))))))
seqtab   <- seqtab[, abs(nchar(colnames(seqtab)) - len_mode) <= len_tolerance, drop = FALSE]

seqtab_nochim <- removeBimeraDenovo(seqtab, method = "consensus", multithread = n_threads)

# Reads retained at each step (Table S19)
track <- tibble(sample = samples, run = unname(run_id),
                filtered = filt[samples, "reads.out"]) %>%
  left_join(bind_rows(counts), by = "sample") %>%
  mutate(nonchim = as.numeric(rowSums(seqtab_nochim)[sample]))
write_csv(track, file.path(out_dir, "read_tracking.csv"))


# ------------------------------------------------------------------------------
# 6. ASV names (by total abundance), FASTA and ASV table
# ------------------------------------------------------------------------------

seqtab_nochim <- seqtab_nochim[, order(colSums(seqtab_nochim), decreasing = TRUE), drop = FALSE]

asv_key <- tibble(asv         = paste0("ASV", seq_len(ncol(seqtab_nochim))),
                  sequence    = colnames(seqtab_nochim),
                  total_reads = colSums(seqtab_nochim))
write_csv(asv_key, file.path(out_dir, "asv_key.csv"))

fa <- DNAStringSet(asv_key$sequence)
names(fa) <- asv_key$asv
writeXStringSet(fa, file.path(out_dir, "asvs.fasta"))

colnames(seqtab_nochim) <- asv_key$asv
sheet <- read_csv(sample_sheet, show_col_types = FALSE) %>%
  filter(seqtab_name %in% rownames(seqtab_nochim))
asv_table <- seqtab_nochim[sheet$seqtab_name, , drop = FALSE]
rownames(asv_table) <- sheet$Sample_ID

# Long format (one row per sample x ASV), used in 004_amf_tables.R
to_long <- function(m) {
  i <- which(m > 0, arr.ind = TRUE)
  tibble(Sample_ID = rownames(m)[i[, 1]], asv = colnames(m)[i[, 2]], reads = m[i]) %>%
    left_join(sheet %>% select(Sample_ID, Species, Sample_type, Batch), by = "Sample_ID")
}
write_csv(to_long(asv_table), file.path(out_dir, "asv_long_nonrarefied.csv"))


# ------------------------------------------------------------------------------
# 7. Rarefaction of all ASVs, per sample type and dataset
# ------------------------------------------------------------------------------

datasets   <- list(batch1 = 1, batch1_2 = c(1, 2))
rare_depth <- list()

for (st in unique(sheet$Sample_type)) for (ds in names(datasets)) {
  ids <- sheet$Sample_ID[sheet$Sample_type == st & sheet$Batch %in% datasets[[ds]]]
  m   <- asv_table[ids, , drop = FALSE]
  depth <- min(rowSums(m))
  set.seed(rare_seed)
  r <- rrarefy(m, sample = depth)
  r <- r[, colSums(r) > 0, drop = FALSE]
  write.csv(t(r), file.path(out_dir, "rarefied", paste0("asv_table_rarefied_", st, "_", ds, ".csv")))
  rare_depth[[paste(st, ds)]] <- tibble(sample_type = st, dataset = ds,
                                        n_samples = nrow(r), rarefied_to = depth)
}

write_csv(bind_rows(rare_depth), file.path(out_dir, "rarefaction_depths.csv"))
