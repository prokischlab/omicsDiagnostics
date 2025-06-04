# eddids: Proteomics analytical & clinical validation

source("src/config.R")
library(matrixStats)
library(stringr)

# 1. Load data ------------------------------------------------------------
sa_file <- "../omicsDagnostics_data/raw_data/proteomics_annotation.tsv"
pg_file <- "../omicsDagnostics_data/raw_data/proteinGroups.txt"

sa <- fread(sa_file)
pg <- fread(pg_file)

# keep bridge channel information
bridge_ids <- sa[NORMALIZATION_SAMPLE == TRUE, SAMPLE_ID]

# Clean proteinGroups similar to MaxQuant defaults
pg <- pg[Reverse != "+" & `Potential contaminant` != "+" & `Only identified by site` != "+"]
pg[, ProteinID := tstrsplit(`Protein IDs`, ";", fixed = TRUE)[[1]]]
setkey(pg, ProteinID)

# 2. Extract reporter intensity and peptide count matrices ----------------
rep_cols <- grep("^Reporter intensity corrected", names(pg), value = TRUE)
pep_cols <- grep("^Razor \+ unique peptides", names(pg), value = TRUE)

rep_mat <- as.matrix(pg[, ..rep_cols])
colnames(rep_mat) <- str_extract(rep_cols, "(?<= )[^ ]+$")
pep_mat <- as.matrix(pg[, ..pep_cols])
colnames(pep_mat) <- str_extract(pep_cols, "(?<= )[^ ]+$")

# 3. Helper to compute metrics -------------------------------------------
get_stats <- function(mat, idx) {
  ok <- !is.na(mat[, idx, drop = FALSE]) & mat[, idx, drop = FALSE] > 0
  data.table(
    Mean   = rowMeans(mat[, idx, drop = FALSE], na.rm = TRUE),
    CV     = rowSds(mat[, idx, drop = FALSE], na.rm = TRUE) /
              rowMeans(mat[, idx, drop = FALSE], na.rm = TRUE),
    Detect = rowMeans(ok),
    Pept   = rowMeans(pep_mat[, idx, drop = FALSE], na.rm = TRUE)
  )
}

stats_bridge <- get_stats(rep_mat, bridge_ids)
stats_all    <- get_stats(rep_mat, colnames(rep_mat))

qc_tbl <- data.table(
  ProteinID       = pg$ProteinID,
  Mean_Int_Bridge = stats_bridge$Mean,
  CV_Bridge       = stats_bridge$CV,
  Pept_Bridge     = stats_bridge$Pept,
  DetRate_Bridge  = stats_bridge$Detect,
  Mean_Int_All    = stats_all$Mean,
  CV_All          = stats_all$CV,
  Pept_All        = stats_all$Pept,
  DetRate_All     = stats_all$Detect
)

qc_tbl[, Tier := fifelse(
  DetRate_Bridge >= 0.9 & CV_Bridge < 0.25 & Pept_Bridge >= 2, "Tier-1",
  fifelse(DetRate_Bridge >= 0.5 & CV_Bridge < 0.40 & Pept_Bridge >= 2, "Tier-2",
          "Tier-3"))]

fwrite(qc_tbl, "QC_metrics_proteins.tsv", sep = "\t")

# 4. Build bridge-based reference sets for analytical validation ----------
ppm <- sweep(rep_mat[, bridge_ids, drop = FALSE], 2,
             colSums(rep_mat[, bridge_ids, drop = FALSE], na.rm = TRUE), "/") * 1e6

prot_detect <- ppm >= 4 & !is.na(ppm)
bridge_pairs <- split(bridge_ids, sa[match(bridge_ids, SAMPLE_ID), PROTEOMICS_BATCH])

det_in_pair <- sapply(bridge_pairs, function(ch) rowSums(prot_detect[, ch, drop = FALSE]) == 2)
num_batches <- rowSums(det_in_pair)

positive <- num_batches >= ceiling(0.8 * length(bridge_pairs))
negative <- num_batches == 0
truth_tbl <- data.table(ProteinID = pg$ProteinID,
                        Truth = fifelse(positive, 1L, fifelse(negative, 0L, NA_integer_)))[!is.na(Truth)]

p_detect <- data.table(ProteinID = pg$ProteinID,
                       p_detect = rowMeans(prot_detect, na.rm = TRUE))
bench <- merge(p_detect, truth_tbl, by = "ProteinID")

# 5. Sensitivity/Specificity across thresholds ---------------------------
confusion_stats <- function(truth, pred) {
  TP <- sum(truth == 1 & pred == 1)
  FP <- sum(truth == 0 & pred == 1)
  FN <- sum(truth == 1 & pred == 0)
  TN <- sum(truth == 0 & pred == 0)
  list(Sens = TP/(TP+FN), Spec = TN/(TN+FP))
}

thr_grid <- seq(0.1, 0.9, by = 0.1)
perf_tbl <- rbindlist(lapply(thr_grid, function(th) {
  pred <- as.integer(bench$p_detect >= th)
  as.data.table(confusion_stats(bench$Truth, pred))[, Threshold := th]
}))

fwrite(perf_tbl, "analytical_performance.tsv", sep = "\t")

# 6. Paper-ready figures -------------------------------------------------
library(ggplot2)
library(ggbeeswarm)

dir.create("Scripts/MultiOmics/Proteomics/clinical_validation/figures", showWarnings = FALSE)

# Scatter of bridge CV vs. mean intensity
p_cv <- ggplot(qc_tbl, aes(x = log10(Mean_Int_Bridge + 1), y = CV_Bridge, colour = Tier)) +
  geom_point(alpha = 0.6) +
  geom_hline(yintercept = 0.25, linetype = "dashed") +
  scale_color_ptol() +
  theme_classic() +
  labs(x = "log10 mean bridge intensity", y = "CV (bridge)",
       colour = "Tier", title = "Technical robustness of protein quantification")

pdf("Scripts/MultiOmics/Proteomics/clinical_validation/figures/CV_vs_intensity.pdf", width = 6, height = 5, useDingbats = FALSE)
print(p_cv)
dev.off()

# Histogram of detection rate across all samples
p_det <- ggplot(qc_tbl, aes(x = DetRate_All, fill = Tier)) +
  geom_histogram(binwidth = 0.05, position = "identity", alpha = 0.6, colour = "black") +
  theme_classic() +
  scale_fill_ptol() +
  labs(x = "Detection rate (all samples)", y = "Number of proteins", fill = "Tier",
       title = "Detection rate distribution")

pdf("Scripts/MultiOmics/Proteomics/clinical_validation/figures/detection_rate_hist.pdf", width = 6, height = 4, useDingbats = FALSE)
print(p_det)
dev.off()
