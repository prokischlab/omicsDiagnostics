#'---
#' title: Aberrant protein expression with covariate regression + zscore
#' author: Ines Scheller
#' wb:
#'  input:
#'  - protrider_annotation: '`sm config["PROC_DATA"] + "/protrider/protrider_annotation.tsv"`'
#'  - protrider_data: '`sm config["PROC_DATA"] + "/protrider/protrider_data.tsv"`'
#'  output:
#'  - covReg_object: '`sm config["PROC_DATA"] + "/covReg/covReg_obj.rds"`'
#'  - covReg_results: '`sm config["PROC_DATA"] + "/covReg/covReg_results.rds"`'
#'  type: script
#'---


source("src/config.R")
library(keras)
use_condaenv('omicsDiagnostics', required = T)   # conda environment with tensorflow 1.14.0 and keras 1.1.0

### load all scripts
scripts_needed = c("data_input.R", "data_handling.R")
for(s in scripts_needed) {
    source(file.path("src/functions/PROTRIDER", s) )
}

# Read annotation
sa <- fread(snakemake@input$protrider_annotation) %>% as.data.frame()
rownames(sa) <- sa$SAMPLE_ID
sa$SAMPLE_ID <- NULL


# Read protrider prepared data
#  genes x samples, dataframe
prot <- fread(snakemake@input$protrider_data) %>% as.data.frame()
rownames(prot) <- prot$geneID
prot$geneID <- NULL

############################################
### read data and annotation as summarizedExperiment object
confounders_used <- c("gender", "PROTEOMICS_BATCH","BATCH_RUN", "INSTRUMENT")
se = summarize_prot(
    file_proteins = prot,
    file_annotation = sa,
    confounders_used = confounders_used )


### regress out covariates
# log2intens <- t(assay(se, "X"))
log2intens <- t(assay(se, "X") + assay(se, "X_mean"))
log2intens[t(assay(se, "X_na"))] <- NA

confounders <- sa[, confounders_used]
confounders$gender[confounders$gender == "male"] <- 0
confounders$gender[confounders$gender == "female"] <- 1
confounders <- as.matrix(confounders)
mode(confounders) <- "numeric"

# covReg <- lm(as.matrix(log2intens) ~ confounders, na.action=na.exclude) # fails because of NAs
# covReg
# resid <- residuals(covReg)
# resid[t(assay(se, "X_na"))] <- NA

models <- apply(log2intens, 2, function(x){
    return(lm(x ~ confounders, na.action=na.exclude))
})
pred <- sapply(models, predict)
resid <- sapply(models, residuals)
zscores <- scale(resid, center=TRUE, scale=TRUE)
summary(rowSums(abs(zscores) > 3, na.rm=TRUE))

covReg_model <- list(models = models, input=log2intens,
                    predicted=pred, residuals=resid, zscores=zscores, 
                    X_mean=t(assay(se, "X_mean")))

############################################
### output and save all the stuff
saveRDS(covReg_model, snakemake@output$covReg_object)


### create results table
# columns:  geneID SAMPLE_ID PROTEIN_LOG2INT PROTEIN_ZSCORE PROTEIN_FC PROTEIN_LOG2FC PROTEIN_PVALUE PROTEIN_PADJ PROTEIN_outlier

long_zscore <- reshape2::melt(zscores, id=colnames(zscores), measure.vars=rownames(zscores))  # value.name does not work

outDf = do.call("cbind", list(long_zscore, 
                              as.vector(log2intens)
                              ))
colnames(outDf) = c('SAMPLE_ID','geneID','PROTEIN_ZSCORE','PROTEIN_LOG2INT')
outDf <- as.data.table(outDf)
# outDf[,PROTEIN_PRED_LOG2INT := as.vector(pred)]
outDf[,PROTEIN_FC := NA]
outDf[,PROTEIN_LOG2FC := NA]
outDf[,PROTEIN_PVALUE:=2*pmin(pnorm(PROTEIN_ZSCORE), 1 - pnorm(PROTEIN_ZSCORE))]
outDf[,PROTEIN_PADJ := NA]
outDf[,PROTEIN_outlier := abs(PROTEIN_ZSCORE) > 3]
# outDf[,PROTEIN_outlier_z5 := abs(PROTEIN_ZSCORE) > 5]

outDf$SAMPLE_ID = as.character(outDf$SAMPLE_ID)
outDf$geneID = as.character(outDf$geneID)

outDf <- outDf[!duplicated(outDf), ]
outDf <- outDf[!is.na(outDf$PROTEIN_ZSCORE), ]

# WRITE RESULTS
saveRDS(outDf, snakemake@output$covReg_results)
