source("src/config.R")
library(psych)
library(reshape2)
library(ggridges)
library(ggpubr)
library(ggfortify)
library(DescTools)
library(pROC)
#load functions
#source("src/functions/LIMMA/limma_functions.R")


# READ ANNOTATION
#sa <- fread(snakemake@input$sample_annotation)
#raw_data <- fread(snakemake@input$raw_prot) %>% as.data.frame()


sa <- fread("../omicsDagnostics_data/raw_data/proteomics_annotation.tsv")
sa <- sa[ NORMALIZATION_SAMPLE == T]
sa$INSTRUMENT <- as.character(sa$INSTRUMENT)
sa[ , INSTRUMENT := paste0("I", INSTRUMENT)]
sa$PROTEOMICS_BATCH <- as.character(sa$PROTEOMICS_BATCH)
sa[ , PROTEOMICS_BATCH := paste0("B", PROTEOMICS_BATCH)]


raw_data <- fread("../omicsDagnostics_data/raw_data/proteomics_not_normalized.tsv") %>% as.data.frame()
rownames(raw_data) <- raw_data$geneID
raw_data$geneID <- NULL
raw_data <- raw_data[ , sa$SAMPLE_ID ]
raw_data[raw_data == 0 ] <- NA # keep missing values as 0 


##############################################################################
## 0.  HELPER FUNCTIONS  ---------------------------------------------------
##############################################################################

## Protein detection call from a numeric matrix + threshold
make_detect_mat <- function(mat, thr, log_space = FALSE) {
  if (log_space) mat >= thr else mat >= thr & !is.na(mat)
}

## Sens/spec at chosen threshold
confusion_stats <- function(truth, pred) {
  TP <- sum(truth == 1 & pred == 1)
  FP <- sum(truth == 0 & pred == 1)
  FN <- sum(truth == 1 & pred == 0)
  TN <- sum(truth == 0 & pred == 0)
  list(Sens = TP/(TP+FN), Spec = TN/(TN+FP))
}


##############################################################################
## 1.  Bridge channels  ----------------------------------------------------
##############################################################################

bridge_ids <- sa[NORMALIZATION_SAMPLE == T ]$SAMPLE_ID    
batch_of    <- setNames(as.character(sa$PROTEOMICS_BATCH), sa$SAMPLE_ID)

## sanity: every batch must have TWO bridge channels
# stopifnot(all(table(sa[match(bridge_ids, SAMPLE_ID), PROTEOMICS_BATCH]) == 2))


 


##############################################################################
## 2.  PPM (TPM‐analogue) calculation  -------------------------------------
##############################################################################
bridge_int  <- raw_data[, bridge_ids] # NA where undetected
 
## convert to PPM per channel  (TPM analogue) ---------------
ppm <- sweep(bridge_int, 2, colSums(bridge_int, na.rm = TRUE), "/") * 1e6
 




##############################################################################
## 3.  Construct positive / negative sets  ---------------------------------
##############################################################################


detect_thr_ppm   <- 4          # PPM ≥ 4  ≈ TPM ≥ 1-3
prot_detect_ppm  <- make_detect_mat(ppm, detect_thr_ppm)

# prot_detect_ppm <- !is.na(ppm) & ppm >= detect_thr_ppm   # logical matrix (genes × bridge channels)

# Group the bridge IDs by PROTEOMICS_BATCH
bridge_pairs <- split(bridge_ids, batch_of[bridge_ids])

# Now bridge_pairs is a list of length = # batches, each with exactly 2 SAMPLE_IDs

det_in_pair <- sapply(bridge_pairs, function(ch) {
  # TRUE if a given protein is detected in both channels of that batch
  rowSums(prot_detect_ppm[, ch, drop = FALSE]) == 2
})


# Count in how many batches each protein passes “both‐channels detected”
num_batches <- rowSums(det_in_pair)

# “Positive” if detected in ≥ 80% of batches; “Negative” if in 0 batches
positive <- num_batches >= ceiling(0.80 * length(bridge_pairs))
negative <- num_batches == 0

truth_tbl <- data.table(
  Gene  = rownames(prot_detect_ppm),
  Truth = fifelse(positive, 1L, fifelse(negative, 0L, NA_integer_))
)[!is.na(Truth)]   # drop “gray zone” proteins not unanimously pos/neg

##############################################################################
## 4.  Protein-wise detection probability  ---------------------------------
##############################################################################

detect_prob <- data.table(
  Gene     = rownames(prot_detect_ppm),
  p_detect = rowMeans(prot_detect_ppm, na.rm = TRUE)
)

bench <- merge(detect_prob, truth_tbl, by = "Gene")

##############################################################################
## 5.  ROC + multi-threshold Sens/Spec table  ------------------------------
##############################################################################
# (a) ROC
roc_obj <- pROC::roc(bench$Truth, bench$p_detect, direction = ">")
plot(roc_obj, col = "firebrick", lwd = 2, print.auc = TRUE,
     legacy.axes = TRUE,
     main = "Bridge-sample benchmark ROC (PPM ≥ 4)")
auc_val <- as.numeric(auc(roc_obj))  # store for reporting


# Main ROC plot
 
plot(roc_obj, col = "firebrick", lwd = 2, print.auc = TRUE,
     legacy.axes = TRUE, main = "Bridge benchmark ROC (PPM ≥ 4)")
 

as.numeric(roc_obj$auc)
# or
roc_obj$auc[1] 


roc_posneg <- roc(bench$Truth, bench$p_detect, direction="<")
as.numeric(roc_posneg$auc)  # likely 1.0 
plot(roc_posneg, print.auc=TRUE, col="darkgreen", main="True AUC (direction='<')")


# (b) Sensitivity/Specificity at thresholds 0.1 … 0.9
# Store sensitivity & specificity at 0.1 … 0.9
thr_grid <- seq(0.1, 0.9, by = 0.1)
grid_tbl <- lapply(thr_grid, function(th) {
  pred <- as.integer(bench$p_detect >= th)
  as.data.table(confusion_stats(bench$Truth, pred))[, Threshold := th]
}) |> rbindlist()

# Bar-style Sens/Spec table
ggplot(melt(grid_tbl, id.vars = "Threshold"), aes(factor(Threshold), value,
                                                  fill = variable)) +
  geom_col(position = "dodge") +
  geom_text(aes(label = sprintf("%.2f", value)), vjust = -0.2, size = 3) +
  scale_fill_ptol(labels = c("Sensitivity", "Specificity")) +
  ylim(0, 1) + theme_classic() +
  labs(x = "p_detect threshold", y = "Metric",
       title = "Analytical performance across thresholds (bridge consensus)")
 
# — 2.1 Histogram of p_detect by Truth group ----------------------------

bench_dt <- as.data.table(bench)
bench_dt[, Truth := factor(Truth, levels = c(0,1), labels = c("Negative","Positive"))]

ggplot(bench_dt, aes(x = p_detect, fill = Truth)) +
  geom_histogram(position = "identity", bins = 50, alpha = 0.6) +
  scale_fill_manual(values = c("#F8766D","#00BFC4")) +
  theme_classic() +
  labs(
    x = "p_detect",
    y = "Number of proteins",
    fill = "Bridge truth",
    title = "Distribution of p_detect by True positive vs. negative"
  )


# — 2.2 Density overlay of p_detect -------------------------------------

ggplot(bench_dt, aes(x = p_detect, color = Truth)) +
  geom_density(size = 1) +
  scale_color_manual(values = c("#F8766D","#00BFC4")) +
  theme_classic() +
  labs(
    x = "p_detect",
    y = "Density",
    color = "Bridge truth",
    title = "Density of p_detect (positive vs. negative)"
  )

# — 2.3 Table of Sens/Spec & line plot -----------------------------------

# grid_tbl already has columns: Threshold, Sens, Spec
grid_dt <- copy(grid_tbl)   # from earlier code; ensure it is a data.table

# Melt to long format for ggplot
grid_melt <- melt(
  grid_dt,
  id.vars = "Threshold",
  measure.vars = c("Sens","Spec"),
  variable.name = "Metric",
  value.name = "Value"
)

ggplot(grid_melt, aes(x = Threshold, y = Value, color = Metric)) +
  geom_line(size = 1.2) +
  geom_point(size = 2) +
  scale_color_manual(values = c("Sensitivity" = "#00BFC4", "Specificity" = "#F8766D")) +
  theme_classic() +
  labs(
    x = "p_detect threshold",
    y = "Metric value",
    color = "Metric",
    title = "Sensitivity & Specificity across p_detect thresholds"
  )

# (C)

bench[, .(
  min_p    = min(p_detect),
  median_p = median(p_detect),
  max_p    = max(p_detect)
), by = Truth ]


##############################################################################
## 6.  Point estimate @ threshold 0.30  ------------------------------------
##############################################################################

best_thr <- 0.30
pred30   <- as.integer(bench$p_detect >= best_thr)
stats30  <- confusion_stats(bench$Truth, pred30)

cat(sprintf(
  "AUC = %.3f |  Sens@0.30 = %.3f | Spec@0.30 = %.3f\n",
  as.numeric(auc(roc_obj)), stats30$Sens, stats30$Spec))

# AUC = 0.000 |  Sens@0.30 = 1.000 | Spec@0.30 = 1.000

##############################################################################
## 7.  Batch-wise sens/spec (non-bridge channels)   ------------------------
##############################################################################

batch_samples <- split(sa$SAMPLE_ID, sa$PROTEOMICS_BATCH)

batch_stats <- rbindlist(lapply(names(batch_samples), function(b) {
  # remove the two bridge channels from that batch
  ch <- setdiff(batch_samples[[b]], bridge_pairs[[b]])
  # “Detected” in that batch if raw intensity > 0 in at least one non-bridge channel
  det_b <- rowSums(!is.na(raw_data[, ch, drop = FALSE]) & raw_data[, ch, drop = FALSE] >= 1) > 0
  tmp   <- data.table(Gene = names(det_b),
                      Pred = as.integer(det_b))
  tmp   <- merge(tmp, truth_tbl, by = "Gene")
  TP <- tmp[Truth == 1 & Pred == 1, .N]
  FP <- tmp[Truth == 0 & Pred == 1, .N]
  FN <- tmp[Truth == 1 & Pred == 0, .N]
  TN <- tmp[Truth == 0 & Pred == 0, .N]
  data.table(Batch = b,
             Sens  = TP/(TP+FN),
             Spec  = TN/(TN+FP))
}))

sens_CI <- MeanCI(batch_stats$Sens)
spec_CI <- MeanCI(batch_stats$Spec)






















############################################################################################################################################################################

# Paper-ready visualisations ----------------------------------------------
library(ggplot2)
library(patchwork)

dir.create("Scripts/MultiOmics/Proteomics/clinical_validation/figures", showWarnings = FALSE)

# ROC curve with ggplot
roc_dt <- data.table(FPR = 1 - roc_obj$specificities,
                     TPR = roc_obj$sensitivities)

p_roc <- ggplot(roc_dt, aes(FPR, TPR)) +
  geom_line(color = "firebrick", size = 1.2) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
  theme_classic() +
  labs(x = "False positive rate", y = "True positive rate",
       title = sprintf("Bridge benchmark ROC (AUC = %.3f)", as.numeric(auc_val)))

pdf("Scripts/MultiOmics/Proteomics/clinical_validation/figures/roc_curve.pdf", width = 5, height = 5, useDingbats = FALSE)
print(p_roc)
dev.off()

# Sensitivity and specificity per batch
p_sens <- ggplot(batch_stats, aes(Batch, Sens)) +
  geom_col(fill = "#00BFC4") +
  geom_hline(yintercept = mean(batch_stats$Sens), linetype = "dashed") +
  theme_classic() +
  labs(x = "Batch", y = "Sensitivity", title = "Sensitivity per batch")

p_spec <- ggplot(batch_stats, aes(Batch, Spec)) +
  geom_col(fill = "#F8766D") +
  geom_hline(yintercept = mean(batch_stats$Spec), linetype = "dashed") +
  theme_classic() +
  labs(x = "Batch", y = "Specificity", title = "Specificity per batch")

batch_plot <- p_sens / p_spec

pdf("Scripts/MultiOmics/Proteomics/clinical_validation/figures/batch_performance.pdf", width = 7, height = 5, useDingbats = FALSE)
print(batch_plot)
dev.off()
