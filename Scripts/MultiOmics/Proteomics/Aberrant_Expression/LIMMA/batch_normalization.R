#'---
#' title: Proteomics data normalization and preparation for limma       
#' author: smirnovd, Chen Meng
#' wb:
#'  input:
#'  - sample_annotation: '`sm config["ANNOTATION"]`'
#'  - raw_prot: '`sm config["RAW_Protein"]`'
#'  output:
#'  - norm_data: '`sm config["PROC_DATA"] + "/limma/proteomics_normalized_not_imputed.tsv"`'
#'  - norm_imp_data: '`sm config["PROC_DATA"] + "/limma/proteomics_normalized_imputed.tsv"`'
#'  - limma_data: '`sm config["PROC_DATA"] + "/limma/limma_data.tsv"`'
#' output: 
#'   html_document:
#'    code_folding: hide
#'    code_download: TRUE
#'---

# source("src/config.R")

#load functions
source("src/functions/LIMMA/limma_functions.R")


# READ ANNOTATION
# sa <- fread('/s/project/mitoMultiOmics/multiOMICs_integration/raw_data/proteomics_annotation.tsv')
sa <- fread(snakemake@input$sample_annotation)




#' # step 1 read and investigate data
step1.data <- fread(snakemake@input$raw_prot) %>% as.data.frame()
# step1.data <- fread('/s/project/mitoMultiOmics/multiOMICs_integration/raw_data/proteomics_not_normalized.tsv') %>% as.data.frame()

rownames(step1.data) <- step1.data$geneID
step1.data$geneID <- NULL

#+ fig.width=10, fig.height=10
plots_normalization(step1.data, sa)


# barplots
dfx <- data.table( SAMPLE_ID = rownames(t(step1.data)) , SUM_INT = rowSums(t(step1.data)) )
dfx <- merge(sa[ ,.(SAMPLE_ID, PROTEOMICS_BATCH, BATCH_RUN, INSTRUMENT)], dfx, by = "SAMPLE_ID" )
dfx <- dfx[order(BATCH_RUN)]

#+ fig.width=12, fig.height=10
ggplot(dfx, aes(SAMPLE_ID , SUM_INT))+
  geom_col()+ 
  theme_bw()+ 
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank())+
  facet_wrap(~PROTEOMICS_BATCH , scales = "free_x", nrow =2 )

rm(dfx)

# Remove samples after inspection
rm.samples <- c("failed_OM54862", "failed_OM12102") 


#' # step 2 between run normalization
# step2.data <- rowwise_normalise(step1.data, quantil_cutoff = F, add.jitter = F, shift1 = FALSE)
step2.data <- rowshift(df = step1.data, batch = sa$PROTEOMICS_BATCH, ref = sa$NORMALIZATION_SAMPLE == T)
#+ fig.width=10, fig.height=10
plots_normalization(step2.data, sa, rm.samples)


#' # step 3 between sample normalisation
step3.data <- colwise_normalise(step2.data)
#+ fig.width=10, fig.height=10
plots_normalization(step3.data, sa, rm.samples)

# Save normalized data
df <- cbind(rownames(step3.data), step3.data)
colnames(df)[1] <- "geneID"
# write_tsv(df, '/s/project/mitoMultiOmics/multiOMICs_integration/processed_data/limma/proteomics_normalized_not_imputed.tsv')
write_tsv(df,  snakemake@output$norm_data)




#' # step 4 Imputation 
step4.data <- impute_min(step2.data, sa)


# Save imputed data
# write_tsv(step4.data, '/s/project/mitoMultiOmics/multiOMICs_integration/processed_data/limma/proteomics_normalized_imputed.tsv')
write_tsv(step4.data,  snakemake@output$norm_imp_data)


genes <- step4.data$geneID
step4.data <- step4.data[ , -1]
rownames(step4.data) <- genes
#+ fig.width=10, fig.height=10
plots_normalization(step4.data, sa, rm.samples)





#' # step 5 Prepare samples for LIMMA (paper cases)
step5.data <- step4.data

sa_paper <- sa[USE_FOR_PROTEOMICS_PAPER == T]

#Generate matrix, subset samples
mat_norm_prot_zero <- as.matrix(step5.data)
mat_norm_prot_zero <- mat_norm_prot_zero[,colnames(mat_norm_prot_zero) %in% sa_paper$SAMPLE_ID ]
dim(mat_norm_prot_zero)
#mat_norm_prot_zero[1:5, 1:5]


## Filter data - removes corrupted samples
# use only samples with less than 30% of NA or 0 
sampleZeroFreq <- apply(mat_norm_prot_zero,2,function(x){ sum(x == 0)/length(x) })
#table(round(sampleZeroFreq, 3))
mat_norm_prot_zero <- mat_norm_prot_zero[, sampleZeroFreq < 0.3] 
dim(mat_norm_prot_zero)


## Filter data - removes proteins with too many NAs
# use only genes expressed in 80% of samples
protZeroFreq <- apply(mat_norm_prot_zero,1,function(x){ sum(x == 0)/length(x) })
table(round(protZeroFreq, 3))
mat_norm_prot_zero <- mat_norm_prot_zero[protZeroFreq < 0.86,]#
dim(mat_norm_prot_zero)
norm_prot_zero <- as.data.frame(mat_norm_prot_zero)


step5.data <- cbind(as.data.frame(rownames(norm_prot_zero)), norm_prot_zero)
setnames(step5.data, "rownames(norm_prot_zero)", "geneID")



# Write data
write_tsv(step5.data,  snakemake@output$limma_data)
# write_tsv(step5.data, '/s/project/mitoMultiOmics/multiOMICs_integration/processed_data/limma/limma_data.tsv')

step5.data$geneID <- NULL
#+ fig.width=10, fig.height=10
plots_normalization(step5.data, sa, rm.samples)
