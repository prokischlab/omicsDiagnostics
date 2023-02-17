#'---
#' title: Proteomics data preparation with imputation for protrider      
#' author: Dmitrii Smirnov
#' wb:
#'  input:
#'  - sa: '`sm config["ANNOTATION"]`'
#'  - raw_prot: '`sm config["RAW_Protein"]`'
#'  output:
#'  - proteomics_imputed: '`sm config["PROC_DATA"] + "/protrider/proteomics_imputed.tsv"`'
#'  - protrider_annotation: '`sm config["PROC_DATA"] + "/protrider/protrider_annotation.tsv"`'
#'  - protrider_data: '`sm config["PROC_DATA"] + "/protrider/protrider_data.tsv"`'
#'  - protrider_data_norm: '`sm config["PROC_DATA"] + "/protrider/protrider_data_norm_log.tsv"`'
#'  type: script
#'---

#load functions
source("src/config.R")
source("src/functions/LIMMA/limma_functions.R")


# READ ANNOTATION
sa <- fread(snakemake@input$sa)



###################
### Impute data ###
###################

#' READ Proteomics data
raw_prot <- fread(snakemake@input$raw_prot) %>% as.data.frame()

rownames(raw_prot) <- raw_prot$geneID
raw_prot$geneID <- NULL

#' Total number of samples
ncol(raw_prot)

# Impute
prot <- impute_min(raw_prot, sa)


# Remove TASP1 gene. This gene has a bad detection rate and was not detected in this sample. 
# Intensity value was imputed.
prot[ prot$geneID == "TASP1",  "OM30476"] <- 0


# Save imputed data
write_tsv(prot,  snakemake@output$proteomics_imputed)



#######################
### Prepare dataset ###
#######################


# Subset samples
sa <- sa[USE_FOR_PROTEOMICS_PAPER ==T]

#' Number of samples
nrow(sa)

# Save subset of the annotation
write_tsv(sa,  snakemake@output$protrider_annotation)


# Subset samples
rownames(prot) <- prot$geneID
prot$geneID <- NULL

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
# table(round(sampleZeroFreq, 3))
mat_prot_zero <- mat_prot_zero[, sampleZeroFreq < 0.3] 
dim(mat_prot_zero)
# heatmap.2(cor(log2(mat_prot_zero + 2)), trace = "none")


#' Filter data - removes proteins with too many NAs
# use only genes expressed in 80% of samples
protZeroFreq <- apply(mat_prot_zero,1,function(x){ sum(x == 0)/length(x) })
# table(round(protZeroFreq, 3))

# MIN_SAMPLE_NUM_PROT <- 0.52 #0.76

print(paste0('proteins detected in >', round(ncol(mat_prot_zero) * (1- MIN_SAMPLE_NUM_PROT)),' samples (', (1- MIN_SAMPLE_NUM_PROT),'%):') )
mat_prot_zero <- mat_prot_zero[protZeroFreq < MIN_SAMPLE_NUM_PROT,]#
dim(mat_prot_zero)
# heatmap.2(cor(log2(mat_prot_zero + 2)), trace = "none")


# Transpose
protrider_data <- cbind( rownames(mat_prot_zero), as.data.frame(mat_prot_zero))
colnames(protrider_data)[1] <- "geneID"
# protrider_data[1:5, 1:5]


# Save data ready for outrider2
write_tsv(protrider_data,  snakemake@output$protrider_data)



# For py_outrider / outrider2

# Set missing values as NA
mat_prot_zero[mat_prot_zero == 0] <- NaN

# Size factor normalisation
sf= DESeq2::estimateSizeFactorsForMatrix(mat_prot_zero)
#barplot(sf, las=2)
# normalize data
mat_norm_prot  <- t(as.matrix(mat_prot_zero))/sf

# Log transformation
norm_prot_log <- as.data.frame( log(mat_norm_prot) )


# Transpose
protrider_data_norm <- cbind( rownames(norm_prot_log), norm_prot_log)
colnames(protrider_data_norm)[1] <- "SAMPLE_ID"
# protrider_data_norm[1:5, 1:5]

# Save norm data
write_tsv(protrider_data_norm,  snakemake@output$protrider_data_norm)

