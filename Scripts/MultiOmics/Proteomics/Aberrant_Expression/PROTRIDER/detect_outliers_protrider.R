#'---
#' title: Aberrant protein expression with PROTRIDER 
#' author: Dmitrii Smirnov
#' wb:
#'  threads: 8
#'  input:
#'   - protrider_config: '`sm config["PROC_DATA"] + "/protrider/config_protrider.yaml"`'
#'   - protrider_annotation: '`sm config["PROC_DATA"] + "/protrider/protrider_annotation.tsv"`'
#'   - protrider_data: '`sm config["PROC_DATA"] + "/protrider/protrider_data.tsv"`'
#'  output:
#'   - protrider_results: '`sm config["PROC_DATA"] + "/protrider/PROTRIDER_results.rds"`'
#'  type: script
#'---

library(data.table)
library(reticulate)

# Use the conda environment defined in wbuild.yaml (via CONDA_ENV) #'   - env_conda: '`sm config["CONDA_ENV"]`'
# use_condaenv(snakemake@input$env_conda, required = TRUE)
use_condaenv('omicsDiagnosticsMinimal', required = T)   

# Define file paths provided by Snakemake
config_path <- snakemake@input$protrider_config
input_intensities <- snakemake@input$protrider_data
sample_annotation <- snakemake@input$protrider_annotation

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
