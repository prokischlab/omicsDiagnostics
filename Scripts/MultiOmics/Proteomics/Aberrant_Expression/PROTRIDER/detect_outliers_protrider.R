#'---
#' title: Aberrant protein expression with PROTRIDER
#' author: loipfins, smirnovd
#' wb:
#'  input:
#'  - protrider_annotation: '`sm config["PROC_DATA"] + "/protrider/protrider_annotation.tsv"`'
#'  - protrider_data: '`sm config["PROC_DATA"] + "/protrider/protrider_data.tsv"`'
#'  output:
#'  - protrider_object: '`sm config["PROC_DATA"] + "/protrider/protrider_obj.rds"`'
#'  - protrider_results: '`sm config["PROC_DATA"] + "/protrider/PROTRIDER_results.rds"`'
#'  type: script
#'---


source("src/config.R")
library(reticulate)
library(keras)
library(dplyr)

# use_condaenv('genetic_diagnosis4', required = T)   # conda environment with tensorflow 1.14.0 and keras 1.1.0
use_condaenv('omicsDiagnostics', required = T)   # conda environment with tensorflow 1.14.0 and keras 1.1.0

### load all scripts
scripts_needed = c("config_protrider.R", "data_input.R", "data_handling.R", "methods_noise_outlier.R", "methods_statistics.R", "autoencoder_models.R")
for(s in scripts_needed) {
  source(file.path("src/functions/PROTRIDER", s) )
}




# Read annotation
# sa <- fread('/s/project/mitoMultiOmics/multiOMICs_integration/processed_data/protrider/protrider_annotation.tsv') %>% as.data.frame()
sa <- fread(snakemake@input$protrider_annotation) %>% as.data.frame()
rownames(sa) <- sa$SAMPLE_ID
sa$SAMPLE_ID <- NULL


# Read protrider prepared data
#  genes x samples, dataframe
prot <- fread(snakemake@input$protrider_data) %>% as.data.frame()
# prot <- fread('/s/project/mitoMultiOmics/multiOMICs_integration/processed_data/protrider/protrider_data.tsv') %>% as.data.frame()
rownames(prot) <- prot$geneID
prot$geneID <- NULL



############################################
### read data and annotation as summarizedExperiment object
se = summarize_prot(
  file_proteins = prot,
  file_annotation = sa,
  confounders_used = c("gender", "PROTEOMICS_BATCH","BATCH_RUN", "INSTRUMENT") )



### define specified metadata in se object
metadata(se)$encod_dim = ENCODING_DIM
metadata(se)$epochs = EPOCHS
metadata(se)$noise = NOISE_FACTOR
metadata(se)$seed = SEED


set.seed(SEED)

############################################
### noise insertion w/o artificial outliers
se = inject_noise(se, NOISE_FACTOR, show_plot = F)


############################################
### autoencoder with predicted values
ae_obj = run_autoencoder_model_wrapper(
  "AE_batches_Dfit",
  se,
  encoding_dim = ENCODING_DIM,
  seed = SEED,
  epochs = EPOCHS,
  batch_size = BATCH_SIZE,
  es_patience = EARLY_STOPPING_PATIENCE,
  optimizer = OPTIMIZER,
  verbose = AE_VERBOSE,
  learning_rate = LEARNING_RATE
)
se = ae_obj$se


############################################
### get signficant outliers
se = calculate_pvalues(se)
se = call_sig_outlier(se, fold_change = 0, pval_adj = 0.1)



############################################
### output and save all the stuff
saveRDS(se, snakemake@output$protrider_object)
#saveRDS(se, "/s/project/mitoMultiOmics/multiOMICs_integration/processed_data/protrider/protrider_obj.rds")




### add extra columns
protr = get_prot_sample_list(se)
setnames(protr,
         c( "protein_id", "sample_id", "intensity_norm_log2", "z_score", "fc", "log2fc", 
             "pvalue", "adj_pvalue", "is_outlier" ),
         c( "geneID", "SAMPLE_ID", "PROTEIN_LOG2INT","PROTEIN_ZSCORE", "PROTEIN_FC","PROTEIN_LOG2FC",
            "PROTEIN_PVALUE" , "PROTEIN_PADJ", "PROTEIN_outlier" ))

protr <-  protr[, c( "geneID", "SAMPLE_ID", "PROTEIN_LOG2INT","PROTEIN_ZSCORE", "PROTEIN_FC","PROTEIN_LOG2FC",
             "PROTEIN_PVALUE" , "PROTEIN_PADJ", "PROTEIN_outlier" )]

protr <- protr[!duplicated(protr), ]
protr <- protr[!is.na(protr$PROTEIN_PVALUE), ]

# WRITE RESULTS
saveRDS(protr, snakemake@output$protrider_results)
# saveRDS(protr,  '/s/project/mitoMultiOmics/multiOMICs_integration/processed_data/protrider/PROTRIDER_results.rds')
