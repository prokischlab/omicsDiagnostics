#'---
#' title: Proteomics data preparation with imputation for PROTRIDER      
#' author: Dmitrii Smirnov
#' wb:
#'  input:
#'   - sa: '`sm config["ANNOTATION"]`'
#'   - raw_prot: '`sm config["RAW_Protein"]`'
#'   - config_protrider: '`sm config["PROTRIDER_CONFIG"]`'
#'   - protr_out: '`sm config["PROC_DATA"] + "/protrider"`'
#'  output:
#'   - proteomics_imputed: '`sm config["PROC_DATA"] + "/protrider/proteomics_imputed.tsv"`'
#'   - protrider_annotation: '`sm config["PROC_DATA"] + "/protrider/protrider_annotation.tsv"`'
#'   - protrider_data: '`sm config["PROC_DATA"] + "/protrider/protrider_data.tsv"`'
#'   - protrider_config: '`sm config["PROC_DATA"] + "/protrider/config_protrider.yaml"`'
#'  type: script
#'---

library(data.table)
library(readr)
library(dplyr)
library(yaml)

# Load helper functions if needed
source("src/config.R")
source("src/functions/LIMMA/limma_functions.R")

# Read sample annotation
sa <- fread(snakemake@input$sa)


# Read raw proteomics data
raw_prot <- fread(snakemake@input$raw_prot) %>% as.data.frame()
rownames(raw_prot) <- raw_prot$geneID
raw_prot$geneID <- NULL

# Impute missing values (using your impute_min function)
prot <- impute_min(raw_prot, sa)

# (Optional) Remove problematic gene (e.g., TASP1)
prot[prot$geneID == "TASP1", "OM30476"] <- 0

# Save imputed data
write_tsv(prot, snakemake@output$proteomics_imputed)

# Subset sample annotation to those used for proteomics and save
sa <- sa[USE_FOR_PROTEOMICS_PAPER == TRUE]
write_tsv(sa, snakemake@output$protrider_annotation)

# Prepare the proteomics data matrix for PROTRIDER
# (Note: We perform filtering but do NOT normalize or log-transform since PROTRIDER handles that)
rownames(prot) <- prot$geneID
prot$geneID <- NULL
prot <- prot[, sa$SAMPLE_ID]
prot <- prot[!duplicated(prot), ]

mat_prot <- as.matrix(prot)
# Filter out samples with >30% zeros
sampleZeroFreq <- apply(mat_prot, 2, function(x) sum(x == 0) / length(x))
mat_prot <- mat_prot[, sampleZeroFreq < 0.3]

# Filter out proteins not expressed in enough samples (e.g., proteins missing in > (MIN_SAMPLE_NUM_PROT*100)% samples)
MIN_SAMPLE_NUM_PROT <- 0.52
message('Filtering proteins detected in >', round(ncol(mat_prot) * (1 - MIN_SAMPLE_NUM_PROT)), ' samples')
protZeroFreq <- apply(mat_prot, 1, function(x) sum(x == 0) / length(x))
mat_prot <- mat_prot[protZeroFreq < MIN_SAMPLE_NUM_PROT, ]

# Save the filtered raw intensities as the PROTRIDER input file
protrider_data <- cbind(rownames(mat_prot), as.data.frame(mat_prot))
colnames(protrider_data)[1] <- "geneID"
write_tsv(protrider_data, snakemake@output$protrider_data)


# ---- Update dataset-specific PROTRIDER configuration ----

# Read the repository config file (e.g. config_protrider.yml)
config_in <- snakemake@input$config_protrider
config_list <- yaml::read_yaml(config_in)

# Update the out_dir to point to the processed data directory for PROTRIDER results.
config_list$out_dir <- snakemake@input$protr_out
# Create the output directory if it doesn't exist.
dir.create(snakemake@input$protr_out, recursive = TRUE, showWarnings = FALSE)

# Write the updated configuration file to output
yaml::write_yaml(config_list, snakemake@output$protrider_config)
