#'---
#' title: Aberrant protein expression with PROTRIDER (OUTRIDER2 implementation)
#' author: Stefan Loipfinger, Ines Scheller, Dmitrii Smirnov
#' wb:
#'  threads: 40
#'  input:
#'  - protrider_object: '`sm config["PROC_DATA"] + "/protrider/protrider_outrider2_obj_unfitted.rds"`'
#'  output:
#'  - protrider_results: '`sm config["PROC_DATA"] + "/protrider/PROTRIDER_results.rds"`'
#'  - protrider_object: '`sm config["PROC_DATA"] + "/protrider/protrider_outrider2_obj.rds"`'
#'  type: script
#'---


# source("src/config.R")
library(reticulate)
use_condaenv('omicsDiagnostics', required = T)   

# load OUTRIDER2 (currently in outrider2 branch of OUTRIDER package)
devtools::load_all("OUTRIDER/") 
register(MulticoreParam(workers=as.integer(snakemake@threads)))

# load Outrider2DataSet
ods <- readRDS(snakemake@input$protrider_object)

## PROTRIDER fit
# confounders that are explicitly considered
confounders_used <- c("gender", "PROTEOMICS_BATCH","BATCH_RUN", "INSTRUMENT")
enc_dim <- getBestQ(ods)
preproParams <- getDefaultPreproParams(ods)
preproParams$noise_factor <- 0.5
ods <- OUTRIDER(ods, q=enc_dim, covariates=confounders_used, seed=111,
                  usePython=TRUE, useBasilisk=FALSE,
                  prepro_options=preproParams, BPPARAM = SerialParam())
ods

############################################
### output and save all the stuff
saveRDS(ods, snakemake@output$protrider_object)



### extract results
protr_res <- results(ods, padjCutoff=0.1, l2fcCutoff=NULL, all=TRUE)
protr_res[, log2fc := preprocessed_raw - preprocessed_expected] # preprocessed values are log2(intensity)
protr_res[, fc := 2^(log2fc)]

# add extra columns
setnames(protr_res,
         c( "featureID", "sampleID", "preprocessed_raw", "preprocessed_expected", "normalized", "zScore", "fc", "log2fc",
            "pValue", "padjust", "aberrant" ),
         c( "geneID", "SAMPLE_ID", "PROTEIN_LOG2INT_RAW", "PROTEIN_LOG2INT_PRED", "PROTEIN_LOG2INT","PROTEIN_ZSCORE", "PROTEIN_FC","PROTEIN_LOG2FC",
            "PROTEIN_PVALUE" , "PROTEIN_PADJ", "PROTEIN_outlier" ))

protr_res <-  protr_res[, c( "geneID", "SAMPLE_ID", "PROTEIN_LOG2INT", # "PROTEIN_LOG2INT_RAW", "PROTEIN_LOG2INT_PRED", 
                             "PROTEIN_ZSCORE", "PROTEIN_FC","PROTEIN_LOG2FC",
                             "PROTEIN_PVALUE" , "PROTEIN_PADJ", "PROTEIN_outlier" )]

protr_res <- protr_res[!duplicated(protr_res), ]
protr_res <- protr_res[!is.na(protr_res$PROTEIN_PVALUE), ]

# WRITE RESULTS
saveRDS(protr_res, snakemake@output$protrider_results)

