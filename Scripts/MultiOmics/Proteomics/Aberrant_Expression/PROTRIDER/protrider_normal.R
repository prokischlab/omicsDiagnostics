#'---
#' title: Protrider Normal  
#' author: Dmitrii Smirnov
#' wb:
#'  input:
#'   - sa: '`sm config["ANNOTATION"]`'
#'   - raw_prot: '`sm config["RAW_DATA"] + "/proteomics_not_normalized_absolute.tsv"`'
#'   - config_protrider: '`sm config["PROTRIDER_CONFIG"]`'
#'   - protr_out: '`sm config["PROC_DATA"] + "/protrider"`'
#'  output:
#'   - protrider_annotation: '`sm config["PROC_DATA"] + "/protrider/protrider_annotation_gaus.tsv"`'
#'   - protrider_data: '`sm config["PROC_DATA"] + "/protrider/protrider_data_gaus.tsv"`'
#'   - protrider_config: '`sm config["PROC_DATA"] + "/protrider/config_protrider_gaus.yaml"`'
#'   - protrider_results: '`sm config["PROC_DATA"] + "/protrider/PROTRIDER_results_gaus.rds"`'
#'  type: script
#'---

library(data.table)
library(readr)
library(dplyr)
library(yaml)




# Load helper functions if needed
source("src/config.R")
source("src/functions/LIMMA/limma_functions.R")

# Create the output directory if it doesn't exist.
dir.create(snakemake@input$protr_out, recursive = TRUE, showWarnings = FALSE)
#dir.create("/Users/Mitya/Desktop/working/omicsDagnostics_data/processed_data/protrider", recursive = TRUE, showWarnings = FALSE)


# Read sample annotation
sa <- fread(snakemake@input$sa)
#sa <- fread("/Users/Mitya/Desktop/working/omicsDagnostics_data/raw_data/proteomics_annotation.tsv")

# Read raw proteomics data
raw_prot <- fread(snakemake@input$raw_prot) %>% as.data.frame()
#raw_prot <- fread("/Users/Mitya/Desktop/working/omicsDagnostics_data/raw_data/proteomics_not_normalized0.tsv") %>% as.data.frame()
rownames(raw_prot) <- raw_prot$geneID
raw_prot$geneID <- NULL
raw_prot[1:5, 1:5]


saZ <- fread("/Users/Mitya/Desktop/working/omicsDagnostics_data/extra/rstudio-export/proteomics_annotation_QC.tsv")
saZ <- saZ[PROTEOME_ID %in% colnames(raw_prot) ]
saZ <- saZ[SAMPLE_ID %in% sa$SAMPLE_ID ]
saZ <- saZ[PROTEOME_ID %in% colnames(raw_prot) ]
sa <- sa[ SAMPLE_ID %in% saZ$SAMPLE_ID ]

setnames(raw_prot, saZ$PROTEOME_ID, saZ$SAMPLE_ID  )


sa <- sa[ SAMPLE_ID %in% saZ$SAMPLE_ID]
sa <- sa[ SAMPLE_ID %in% colnames(raw_prot) ]
raw_prot <- raw_prot[ , sa$SAMPLE_ID]


# Impute missing values (using your impute_min function)
prot <- impute_min(raw_prot, sa)

 

 

# Subset sample annotation to those used for proteomics and save
sa <- sa[USE_FOR_PROTEOMICS_PAPER == TRUE]
write_tsv(sa, snakemake@output$protrider_annotation)
# write_tsv(sa, "/Users/Mitya/Desktop/working/omicsDagnostics_data/processed_data/protrider/protrider_annotation_gaus.tsv")

 
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

mat_prot[1:5, 1:5]

# # Filter out proteins not expressed in enough samples (e.g., proteins missing in > (MIN_SAMPLE_NUM_PROT*100)% samples)
# MIN_SAMPLE_NUM_PROT <- 0.52
# message('Filtering proteins detected in >', round(ncol(mat_prot) * (1 - MIN_SAMPLE_NUM_PROT)), ' samples')
# protZeroFreq <- apply(mat_prot, 1, function(x) sum(x == 0) / length(x))
# mat_prot <- mat_prot[protZeroFreq < MIN_SAMPLE_NUM_PROT, ]

# Save the filtered raw intensities as the PROTRIDER input file
protrider_data <- cbind(rownames(mat_prot), as.data.frame(mat_prot))
colnames(protrider_data)[1] <- "geneID"
write_tsv(protrider_data, snakemake@output$protrider_data)
#write_tsv(protrider_data, "/Users/Mitya/Desktop/working/omicsDagnostics_data/processed_data/protrider/protrider_data_gaus.tsv")

# ---- Update dataset-specific PROTRIDER configuration ----

# Read the repository config file (e.g. config_protrider.yml)
config_in <- snakemake@input$config_protrider
config_in <- "config_protrider.yml"
config_list <- yaml::read_yaml(config_in)

# Update the out_dir to point to the processed data directory for PROTRIDER results.
config_list$out_dir <- snakemake@input$protr_out
#config_list$out_dir <- "/Users/Mitya/Desktop/working/omicsDagnostics_data/processed_data/protrider"

config_list$pval_dist <- "gaussian"

# Write the updated configuration file to output
yaml::write_yaml(config_list, snakemake@output$protrider_config)
#yaml::write_yaml(config_list, "/Users/Mitya/Desktop/working/omicsDagnostics_data/processed_data/protrider/config_protrider_gaus.yaml")


library(data.table)
library(reticulate)



# Use the conda environment defined in wbuild.yaml (via CONDA_ENV) #'   - env_conda: '`sm config["CONDA_ENV"]`'
# use_condaenv(snakemake@input$env_conda, required = TRUE)
use_condaenv('omicsDiagnosticsMinimal', required = T)   

# Define file paths provided by Snakemake
config_path <- snakemake@input$protrider_config
input_intensities <- snakemake@output$protrider_data
sample_annotation <- snakemake@output$protrider_annotation



# config_path <- "/Users/Mitya/Desktop/working/omicsDagnostics_data/processed_data/protrider/config_protrider_gaus.yaml"
# input_intensities <- "/Users/Mitya/Desktop/working/omicsDagnostics_data/processed_data/protrider/protrider_data_gaus.tsv"
# sample_annotation <- "/Users/Mitya/Desktop/working/omicsDagnostics_data/processed_data/protrider/protrider_annotation_gaus.tsv"


# Build and run the PROTRIDER command-line call
cmd <- "protrider"
args <- c("--config", config_path,
          "--input_intensities", input_intensities,
          "--sample_annotation", sample_annotation)

message("Running PROTRIDER with command:\n", paste(cmd, args, collapse=" "))
system2(cmd, args = args)

message("Protrider run finished.")

# Load the results (assumed to be written as a CSV summary file in the out_dir specified by config)
config_list <- yaml::read_yaml(config_path)
result_path <- file.path(config_list$out_dir, "/protrider_summary.csv")
res <- fread(result_path)

setnames(res, c("sampleID", "proteinID"), c("SAMPLE_ID", "geneID"))

saveRDS(res, snakemake@output$protrider_results)
 


 

