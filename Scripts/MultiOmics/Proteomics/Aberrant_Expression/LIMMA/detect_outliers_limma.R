#'---
#' title: Aberrant protein expression with LIMMA
#' author: Dmitrii Smirnov
#' wb:
#'  input:
#'  - sa: '`sm config["ANNOTATION"]`'
#'  - limma_data: '`sm config["PROC_DATA"] + "/limma/limma_data.tsv"`'
#'  output:
#'  - limma: '`sm config["PROC_DATA"] + "/limma/LIMMA_results.rds"`'
#'  type: script
#'---


# load functions
#source("src/config.R")
source("src/functions/LIMMA/limma_functions.R")



#Read proteomics data
limma_data <- fread(snakemake@input$limma_data) %>% as.data.frame()

genes <- limma_data$geneID
limma_data <- as.matrix(limma_data[ , -1])
rownames(limma_data) <- genes
dim(limma_data)

#Read sample annotation
sa <- fread(snakemake@input$sa)
sa <- sa[USE_FOR_PROTEOMICS_PAPER == T]

sa[gender == 'male', gender := 1]
sa[gender == 'female', gender := 0]
sa[INSTRUMENT == 1, INSTRUMENT := 0]
sa[INSTRUMENT == 2, INSTRUMENT := 1]


#Subset annotation for the analysis
sa <- sa[, c("SAMPLE_ID", "gender", "PROTEOMICS_BATCH", "INSTRUMENT", "BATCH_RUN") ]
sa <- sa[!duplicated(sa), ]
sa$gender <- as.numeric(sa$gender)
rownames(sa) <- sa$SAMPLE_ID
sa$SAMPLE_ID <- NULL


# run aberrant protein abundances analysis
prot <- wrapper_aberrant_protein_expr_simple(prot_intensity= limma_data, 
                                             coln_sample_id = "SAMPLE_ID", 
                                             p_adjust_method = "BY", 
                                             LIMMA_covariates = sa,
                                             normalize_matrix = T) 

setnames(prot, "GENEID", "geneID")
prot$PROTEIN_INT <- 2^(prot$PROTEIN_LOG2INT)
prot$PROTEIN_FC <- 2^(prot$PROTEIN_LOG2FC)

prot <- prot[, c("geneID" , "SAMPLE_ID", "PROTEIN_LOG2INT","PROTEIN_ZSCORE", "PROTEIN_FC","PROTEIN_LOG2FC",
                 "PROTEIN_PVALUE" , "PROTEIN_PADJ", "PROTEIN_outlier" ) ]


prot[, validated := F]
prot[ PROTEIN_outlier == T &  PROTEIN_LOG2FC < 0, validated := T]
prot[PROTEIN_LOG2FC <= -1 & PROTEIN_PVALUE < 0.05 , validated := T]
prot[PROTEIN_ZSCORE <= -2 & PROTEIN_PVALUE < 0.05  , validated := T]
prot <- prot[!is.na(PROTEIN_PVALUE)]

# WRITE RESULTS
saveRDS(prot,  snakemake@output$limma)











