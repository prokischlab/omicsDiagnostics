#'---
#' title: Aberrant protein expression with PROTRIDER (OUTRIDER2 implementation)
#' author: Stefan Loipfinger, Ines Scheller
#' wb:
#'  threads: 15
#'  input:
#'  - protrider_annotation: '`sm config["PROC_DATA"] + "/protrider/protrider_annotation.tsv"`'
#'  - protrider_data: '`sm config["PROC_DATA"] + "/protrider/protrider_data.tsv"`'
#'  output:
#'  - protrider_object: '`sm config["PROC_DATA"] + "/protrider/protrider_outrider2_obj_unfitted.rds"`'
#'  type: script
#'---


# source("src/config.R")
library(reticulate)
use_condaenv('omicsDiagnostics', required = T)   

# load OUTRIDER2 (currently in outrider2 branch of OUTRIDER package)
devtools::load_all("OUTRIDER/") 
# register(MulticoreParam(workers=4))
register(MulticoreParam(workers=as.integer(snakemake@threads)))

## input files
protrider_annotation <- snakemake@input$protrider_annotation
protrider_data <- snakemake@input$protrider_data

# read input files
sa <- fread(protrider_annotation) 
setnames(sa, "SAMPLE_ID", "sampleID")
sa[gender == "male", gender := 0]
sa[gender == "female", gender := 1]
sa[, gender := as.integer(gender)]

prot <- fread(protrider_data) %>% as.data.frame()
rownames(prot) <- prot$geneID
prot$geneID <- NULL
prot[prot == 0] <- NA

# create Outrider2DataSet
ods <- Outrider2DataSet(inputData=prot, colData=sa, profile="protrider")
preproParams <- getDefaultPreproParams(ods)
preproParams$noise_factor <- 0.5

# ods <- filterExpression(ods, max_na_percentage=0.75) # could be run for filtering genes with too many NAs

## PROTRIDER hyper param opt
# confounders that are explicitly considered
confounders_used <- c("gender", "PROTEOMICS_BATCH","BATCH_RUN", "INSTRUMENT")
ods <- findEncodingDim(ods, covariates=confounders_used, seed=111,
                        usePython=TRUE, useBasilisk=FALSE,
                        prepro_options=preproParams)

plotEncDimSearch(ods)
metadata(ods)$encDimTable


############################################
### output and save all the stuff
saveRDS(ods, snakemake@output$protrider_object)

