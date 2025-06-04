source("src/config.R")
library(psych)
library(reshape2)
library(ggridges)
library(ggpubr)
library(ggfortify)
#load functions
#source("src/functions/LIMMA/limma_functions.R")


# READ ANNOTATION
#sa <- fread(snakemake@input$sample_annotation)
#raw_data <- fread(snakemake@input$raw_prot) %>% as.data.frame()


sa <- fread("../omicsDagnostics_data/raw_data/proteomics_annotation.tsv")

uniqueN(sa$BATCH_RUN)

sa <- sa[ NORMALIZATION_SAMPLE == T]
sa$INSTRUMENT <- as.character(sa$INSTRUMENT)
sa[ , INSTRUMENT := paste0("I", INSTRUMENT)]

sa$PROTEOMICS_BATCH <- as.character(sa$PROTEOMICS_BATCH)
sa[ , PROTEOMICS_BATCH := paste0("B", PROTEOMICS_BATCH)]

sa$BATCH_RUN <- as.character(sa$BATCH_RUN)
sa[ , BATCH_RUN := paste0("BR", BATCH_RUN)]

saX <- as.data.frame(sa[ , c("SAMPLE_ID", "PROTEOMICS_BATCH", "BATCH_RUN", "INSTRUMENT"    )])
rownames(saX) <- saX$SAMPLE_ID
saX$SAMPLE_ID <- NULL
saX$INSTRUMENT <- as.factor(saX$INSTRUMENT)
saX$PROTEOMICS_BATCH <- as.factor(saX$PROTEOMICS_BATCH)
saX$BATCH_RUN <- as.factor(saX$BATCH_RUN)
head(saX)

raw_data <- fread("../omicsDagnostics_data/raw_data/proteomics_not_normalized.tsv") %>% as.data.frame()
rownames(raw_data) <- raw_data$geneID
raw_data$geneID <- NULL

raw_data <- raw_data[ , sa$SAMPLE_ID ]
raw_data[raw_data == 0 ] <- NA

raw_data_not_log <- raw_data
raw_data <- log2(raw_data)
head(raw_data)




raw_dataZ <- raw_data
raw_dataZ[is.na(raw_dataZ)] <- min(raw_data, na.rm = T)
pca_res <- prcomp(t(raw_dataZ), scale. = TRUE, center = TRUE)

autoplot(pca_res, data = saX, colour = 'PROTEOMICS_BATCH', shape = 'INSTRUMENT', frame = TRUE) +
  theme_classic() +
  ggtitle("PCA of Proteomics Data Annotated by Batch and Instrument")






## Pairwise cor 

# Calculate correlations and p-values handling missing values pairwise
corr_result <- corr.test(raw_data, use = "pairwise.complete.obs", method = "pearson", adjust = "none")

pheatmap(corr_result$r, annotation_col  = saX , annotation_row  = saX)


pheatmap(corr_result$r, 
         annotation_col = saX, 
         annotation_row = saX,
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation",
         clustering_method = "ward.D2",
         show_rownames = FALSE,
         show_colnames = FALSE,
         main = "Proteomics Sample Correlation Heatmap")


##############################
### Long Cor
##############################

# No lower/upper.tri manipulation at this stage
cor_mat <- corr_result$r
pval_mat <- corr_result$p

# Convert fully to long format first, including diagonals
cor_long <- melt(cor_mat, varnames = c("Variable1", "Variable2"), value.name = "Correlation")
pval_long <- melt(pval_mat, varnames = c("Variable1", "Variable2"), value.name = "p_value")

# Merge correlations and p-values into one long-format dataframe
result_long <- merge(cor_long, pval_long, by = c("Variable1", "Variable2"))

# Remove self-correlations explicitly
result_long <- result_long[result_long$Variable1 != result_long$Variable2, ]

# Remove duplicates systematically: ensure you only keep unique pairs irrespective of order
result_long$Pair <- apply(result_long[, c("Variable1", "Variable2")], 1, function(x) paste(sort(x), collapse = "_"))
result_long <- result_long[!duplicated(result_long$Pair), ]
result_long$Pair <- NULL  # Clean-up

# Merge annotations clearly and consistently
result_longX <- merge(result_long, sa[, .(SAMPLE_ID, PROTEOMICS_BATCH, BATCH_RUN, INSTRUMENT)],
                      by.x = "Variable1", by.y = "SAMPLE_ID", all.x = TRUE)
result_longX <- merge(result_longX, sa[, .(SAMPLE_ID, PROTEOMICS_BATCH, BATCH_RUN, INSTRUMENT)],
                      by.x = "Variable2", by.y = "SAMPLE_ID", suffixes = c(".x", ".y"), all.x = TRUE)

result_longX <- as.data.table(result_longX)

# Define groups explicitly and clearly:
result_longX[, BATCH := ifelse(PROTEOMICS_BATCH.x == PROTEOMICS_BATCH.y, "same", "different")]
result_longX[, Fractionation_BATCH := ifelse(BATCH_RUN.x == BATCH_RUN.y, "same", "different")]
result_longX[, INSTRUMENT := "different"]
result_longX[INSTRUMENT.x == INSTRUMENT.y & INSTRUMENT.x == "I1", INSTRUMENT := "Instrument1"]
result_longX[INSTRUMENT.x == INSTRUMENT.y & INSTRUMENT.x == "I2", INSTRUMENT := "Instrument2"]

# Clean up final columns
result_longX <- result_longX[ , c("Variable1", "Variable2","Correlation", "p_value",
                                  "BATCH", "Fractionation_BATCH", "INSTRUMENT")]


head(result_longX)















result_longX_summary <- result_longX[, .(
  Median_Corr = median(Correlation, na.rm = TRUE),
  Mean_Corr = mean(Correlation, na.rm = TRUE),
  SD_Corr = sd(Correlation, na.rm = TRUE)
), by = .(BATCH, Fractionation_BATCH, INSTRUMENT)]

print(result_longX_summary)


#View(result_longX[BATCH == "same" ])

# By batch  #########################
result_longX[, .(
  Median_Corr = median(Correlation, na.rm = TRUE),
  Mean_Corr = mean(Correlation, na.rm = TRUE),
  SD_Corr = sd(Correlation, na.rm = TRUE)
), by = .(BATCH)]

ggplot( result_longX, aes(BATCH, Correlation ))+
  geom_quasirandom(aes(color = BATCH))+
  geom_boxplot()+
  stat_compare_means()+
  scale_color_ptol()+
  theme_classic()

ggplot(result_longX, aes(x = Correlation, y = BATCH, fill = BATCH)) +
  geom_density_ridges(alpha = 0.7) +
  theme_classic() +
  scale_fill_ptol() +
  stat_compare_means(label.x = 0.5) +
  ggtitle("Correlation distribution by Batch")



# By fractionation batch  ################  ################  ################
result_longX[, .(
  Median_Corr = median(Correlation, na.rm = TRUE),
  Mean_Corr = mean(Correlation, na.rm = TRUE),
  SD_Corr = sd(Correlation, na.rm = TRUE)
), by = .(  Fractionation_BATCH )]

ggplot( result_longX, aes(Fractionation_BATCH, Correlation ))+
  geom_quasirandom(aes(color = Fractionation_BATCH))+
  geom_boxplot()+
  stat_compare_means()+
  scale_color_ptol()+
  theme_classic()
# Fractionation batch visualization
ggplot(result_longX, aes(x = Correlation, y = Fractionation_BATCH, fill = Fractionation_BATCH)) +
  geom_density_ridges(alpha = 0.7) +
  theme_classic() +
  scale_fill_ptol() +
  stat_compare_means(label.x = 0.5) +
  ggtitle("Correlation distribution by Fractionation Batch")




# By instrument   ################  ################  ################
result_longX[, .(
  Median_Corr = median(Correlation, na.rm = TRUE),
  Mean_Corr = mean(Correlation, na.rm = TRUE),
  SD_Corr = sd(Correlation, na.rm = TRUE)
), by = .( INSTRUMENT)]

# INSTRUMENT Median_Corr Mean_Corr    SD_Corr
# <char>       <num>     <num>      <num>
# 1: Instrument2   0.8846086 0.8898447 0.02933564
# 2:   different   0.8647461 0.8639648 0.01848894
# 3: Instrument1   0.8747273 0.8827541 0.03742202

ggplot( result_longX, aes(INSTRUMENT, Correlation ))+
  geom_boxplot(aes(fill = INSTRUMENT))+
  geom_quasirandom(alpha = 0.5)+
  stat_compare_means()+
  scale_color_ptol()+
  theme_classic()

# Instrument effect visualization
ggplot(result_longX[INSTRUMENT != "different"], aes(x = Correlation, y = INSTRUMENT, fill = INSTRUMENT)) +
  geom_density_ridges(alpha = 0.7) +
  theme_classic() +
  scale_fill_ptol() +
  stat_compare_means(label.x = 0.5) +
  ggtitle("Correlation distribution by Instrument")

# Check result
head(result_longX)

# Todo
# remove duplicates (b to a), keep only a to b
# Annotate with Batches, Instruments etc, Batch run 
# analyse - within batch correlation coefitients 
# analyse - between batch correlation coefitients 
# analyse effect of different instruments _ within instruments, between instruments 







# Batch visualization
ggplot(result_longX, aes(x = Correlation, fill = BATCH)) +
  geom_density(alpha = 0.5) +
  theme_classic() +
  scale_fill_ptol() +
  labs(title = "Density of Correlations (Within vs Between Batches)",
       x = "Correlation", y = "Density")  

# Paper-ready summary figure
library(patchwork)

p_batch <- ggplot(result_longX, aes(x = Correlation, fill = BATCH)) +
  geom_density(alpha = 0.7) +
  scale_fill_ptol() +
  theme_classic() +
  labs(title = "Correlation by Batch", x = "Pearson r", y = "Density")

p_frac <- ggplot(result_longX, aes(x = Correlation, fill = Fractionation_BATCH)) +
  geom_density(alpha = 0.7) +
  scale_fill_ptol() +
  theme_classic() +
  labs(title = "Correlation by Fractionation", x = "Pearson r", y = "Density")

p_instr <- ggplot(result_longX[INSTRUMENT != "different"], aes(x = Correlation, fill = INSTRUMENT)) +
  geom_density(alpha = 0.7) +
  scale_fill_ptol() +
  theme_classic() +
  labs(title = "Correlation by Instrument", x = "Pearson r", y = "Density")

paper_plot <- (p_batch | p_frac) / p_instr

dir.create("Scripts/MultiOmics/Proteomics/clinical_validation/figures", showWarnings = FALSE)

pdf("Scripts/MultiOmics/Proteomics/clinical_validation/figures/reproducibility_overview.pdf", width = 8, height = 6, useDingbats = FALSE)
print(paper_plot)
dev.off()
