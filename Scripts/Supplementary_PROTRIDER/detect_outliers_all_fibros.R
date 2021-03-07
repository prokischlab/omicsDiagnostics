#'---
#' title: Protrider all in one script to detect outliers in the fibroblast samples      
#' author: smirnovd
#' wb:
#'  input:
#'  - sa: '`sm config["ANNOTATION"]`'
#'  - raw_prot: '`sm config["RAW_Protein"]`'
#'  output:
#'  - protrider_results_fib: '`sm config["PROC_DATA"] + "/protrider/results_fib_all.rds"`'
#'  type: script
#'---

#load functions
source("src/config.R")
source("src/functions/LIMMA/limma_functions.R")
library(reticulate)
library(keras)
library(dplyr)

use_condaenv('omicsDiagnostics', required = T)   # conda environment with tensorflow 1.14.0 and keras 1.1.0


### load all scripts
scripts_needed = c("config_protrider.R", "data_input.R", "data_handling.R", "methods_noise_outlier.R", "methods_statistics.R", "autoencoder_models.R")
for(s in scripts_needed) {
  source(file.path("src/functions/PROTRIDER", s) )
}



# READ ANNOTATION
# sa <- fread('/s/project/mitoMultiOmics/multiOMICs_integration/raw_data/proteomics_annotation.tsv')
sa <- fread(snakemake@input$sa)



###################
### Impute data ###
###################

#' READ Proteomics data
raw_prot <- fread(snakemake@input$raw_prot) %>% as.data.frame()
# raw_prot <- fread('/s/project/mitoMultiOmics/multiOMICs_integration/raw_data/proteomics_not_normalized.tsv') %>% as.data.frame()

rownames(raw_prot) <- raw_prot$geneID
raw_prot$geneID <- NULL

#' Total number of samples
ncol(raw_prot)

# Impute
prot <- impute_min(raw_prot, sa)

# Remove TASP1 gene. This gene has a bad detection rate and was not detected in this sample. 
# Intensity value was imputed.
prot[ prot$geneID == "TASP1",  "OM30476"] <- 0

# Subset samples
rownames(prot) <- prot$geneID
prot$geneID <- NULL


#' Don't run imputation 
prot <- raw_prot

rm(raw_prot)

#######################
### Prepare dataset ###
#######################


# Subset samples
sa <- sa[TISSUE == "FIBROBLAST"] #  

# remove duplicated normalisation samples 
sa <- sa[!(NORMALIZATION_SAMPLE == T &  USE_FOR_PROTEOMICS_PAPER == F)] 

# remove replicates
sa <- sa[!( !is.na(REPLICATE) &  USE_FOR_PROTEOMICS_PAPER == F)] 

#' Number of samples
nrow(sa)




# Subset samples
prot <- prot[ , sa$SAMPLE_ID]
prot <- prot[!duplicated(prot ),]


# Create matrix
mat_prot_zero <- as.matrix(prot)
dim(mat_prot_zero)
# mat_prot_zero[1:5, 1:5]


#' Filter data - removes corrupted samples
# use only samples with less than 30% of NA or 0 
sampleZeroFreq <- apply(mat_prot_zero,2,function(x){ sum(x == 0)/length(x) })
table(round(sampleZeroFreq, 3))
mat_prot_zero <- mat_prot_zero[, sampleZeroFreq < 0.3] 
dim(mat_prot_zero)
# heatmap.2(cor(log2(mat_prot_zero + 2)), trace = "none")


#' Filter data - removes proteins with too many NAs
protZeroFreq <- apply(mat_prot_zero, 1, function(x){ sum(x == 0)/length(x) })
table(round(protZeroFreq, 3))

MIN_SAMPLE_NUM_PROT <- 0.50
print(paste0('proteins detected in >', round(ncol(mat_prot_zero) * (1- MIN_SAMPLE_NUM_PROT)),' samples (', (1- MIN_SAMPLE_NUM_PROT),'%):') )
mat_prot_zero <- mat_prot_zero[protZeroFreq < MIN_SAMPLE_NUM_PROT,]#
dim(mat_prot_zero)
# heatmap.2(cor(log2(mat_prot_zero + 2)), trace = "none")




#######################
### Detect outliers ###
#######################
sa <- as.data.frame(sa)
rownames(sa) <- sa$SAMPLE_ID
sa$SAMPLE_ID <- NULL

prot <- as.data.frame(mat_prot_zero)


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
saveRDS(protr, snakemake@output$protrider_results_fib)
# saveRDS(protr,  '/s/project/mitoMultiOmics/multiOMICs_integration/processed_data/protrider/results_fib_all.rds')



