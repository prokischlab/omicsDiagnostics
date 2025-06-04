

#load functions
source("/data/agprokisch/working/dima/proteomics_analysis/src/config.R")
source("/data/agprokisch/working/dima/proteomics_analysis/src/functions/LIMMA/limma_functions.R")

suppressPackageStartupMessages({
  library(reticulate)
  library(ggplot2)
  library(data.table)
  library(dplyr)
  library(magrittr)
  library(readr)
  library(tibble)
  library(SummarizedExperiment)
  library(OUTRIDER)
})


# load OUTRIDER2 (currently in outrider2 branch of OUTRIDER package)



################################################
######## Data prep #############################
################################################


# READ ANNOTATION
sa <- fread('/data/agprokisch/working/data/proteomics/raw_data/proteomics_annotation.tsv')

unique(sa$GROUP)
sa[ GROUP == "", GROUP := NA]
sa[ CATEGORY %in% c("IIb", "IIc"), SOLVED := "NO"]
sa <- sa[USE_FOR_PROTEOMICS_PAPER == T | ( USE == T & GROUP == "RETT_RETT_LIKE" ) ]



res_ens_rett <- readRDS('/data/agprokisch/working/data/proteomics/processed_results/protrider/paper_rett_results_savee.rds') %>% as.data.table()
prot_iterX <- res_ens_rett[ , c("PROTEOME_ID", "geneID", "PROTEIN_LOG2INT")]
saX <- sa[ PROTEOME_ID %in%  unique(prot_iterX$PROTEOME_ID)]

mecp2 <-  prot_iterX[ geneID == "MECP2"]

prot_iterations <- c()


# i <- 70
# for( i in  ( 10 : (uniqueN(saX$PROTEOME_ID) - 1) ) ){
# 
#   for (j in seq( 100 )) {
#     print(paste("N samples:", i, "Iteration:", j ))
#     sa_temp <- saX[sample(nrow(saX), i), ]
# 
#     prot_iter <- mecp2[ PROTEOME_ID %in% sa_temp$PROTEOME_ID ]
#     # prot_iter <- prot_iterX[ PROTEOME_ID %in% sa_temp$PROTEOME_ID ]
#     prot_iter[ , MEAN := fitdistr(na.exclude(PROTEIN_LOG2INT), "normal")$estimate[1] , by = geneID]
#     prot_iter[ , SD := fitdistr(na.exclude(PROTEIN_LOG2INT), "normal")$estimate[2] , by = geneID ]
#     prot_iter[ , PROTEIN_LOG2FC := PROTEIN_LOG2INT - MEAN ]
#     prot_iter[ , PROTEIN_ZSCORE := PROTEIN_LOG2FC/SD ]
#     prot_iter[ , Pval := pnorm(PROTEIN_LOG2INT , mean = MEAN, sd = SD ) ]
#     prot_iter[ , PROTEIN_PVALUE := 2*pmin(Pval, 1-Pval) ]
#     prot_iter$Pval <- NULL
#     prot_iter[ , PROTEIN_PADJ := p.adjust(PROTEIN_PVALUE, method = 'BY'), by = PROTEOME_ID ]
#     prot_iter[ , SampleSize := i ]
#     prot_iter[ , Iteration := j]
#     prot_iterations <- rbind(prot_iterations, prot_iter)
#   }
# }
# 
# 
# saveRDS(prot_iterations,  '/data/agprokisch/working/data/proteomics/processed_results/protrider/iterate_full_dataset_MECP2.rds')


causal <- prot_iterX[ geneID  %in% unique(sa$KNOWN_MUTATION)]
prot_iterations <- c()
# # i <- 70
# for( i in  ( 10 : 220 ) ){
#   
#   for (j in seq( 230 - i ) ) {
#     print(paste("N samples:", i, "Iteration:", j ))
#     sa_temp <- saX[sample(nrow(saX), i), ]
#     
#     prot_iter <- causal[ PROTEOME_ID %in% sa_temp$PROTEOME_ID ]
#     prot_iter[ , MEAN := fitdistr(na.exclude(PROTEIN_LOG2INT), "normal")$estimate[1] , by = geneID]
#     prot_iter[ , SD := fitdistr(na.exclude(PROTEIN_LOG2INT), "normal")$estimate[2] , by = geneID ]
#     prot_iter[ , PROTEIN_LOG2FC := PROTEIN_LOG2INT - MEAN ]
#     prot_iter[ , PROTEIN_ZSCORE := PROTEIN_LOG2FC/SD ]
#     prot_iter[ , Pval := pnorm(PROTEIN_LOG2INT , mean = MEAN, sd = SD ) ] 
#     prot_iter[ , PROTEIN_PVALUE := 2*pmin(Pval, 1-Pval) ] 
#     prot_iter$Pval <- NULL
#     prot_iter[ , PROTEIN_PADJ := p.adjust(PROTEIN_PVALUE, method = 'BY'), by = PROTEOME_ID ]
#     prot_iter[ , SampleSize := i ]
#     prot_iter[ , Iteration := j]
#     prot_iterations <- rbind(prot_iterations, prot_iter)
#   }
# }
# saveRDS(prot_iterations,  '/data/agprokisch/working/data/proteomics/processed_results/protrider/iterate_full_dataset_causal.rds')


saX[is.na(KNOWN_MUTATION), KNOWN_MUTATION := "unsolved" ]
sa_main <- saX[ KNOWN_MUTATION != "MECP2"]
sa_mecp2 <- saX[ KNOWN_MUTATION == "MECP2"]


Nmecp2 <- uniqueN(sa_mecp2$PROTEOME_ID)

causal <- prot_iterX[ geneID  %in% unique(sa$KNOWN_MUTATION)]
prot_iterations <- c()

# i <- 70
for( i in 1: Nmecp2  ){
  
  for (j in seq( 150 )) {
    print(paste("N samples:", i, "Iteration:", j ))

    
    sa_tempM <- sa_mecp2[sample(nrow(sa_mecp2), i), ]
    sa_temp <- rbind(sa_main, sa_tempM)
    
    

    prot_iter <- causal[ PROTEOME_ID %in% sa_temp$PROTEOME_ID ]
    prot_iter[ , MEAN := fitdistr(na.exclude(PROTEIN_LOG2INT), "normal")$estimate[1] , by = geneID]
    prot_iter[ , SD := fitdistr(na.exclude(PROTEIN_LOG2INT), "normal")$estimate[2] , by = geneID ]
    prot_iter[ , PROTEIN_LOG2FC := PROTEIN_LOG2INT - MEAN ]
    prot_iter[ , PROTEIN_ZSCORE := PROTEIN_LOG2FC/SD ]
    prot_iter[ , Pval := pnorm(PROTEIN_LOG2INT , mean = MEAN, sd = SD ) ]
    prot_iter[ , PROTEIN_PVALUE := 2*pmin(Pval, 1-Pval) ]
    prot_iter$Pval <- NULL
    prot_iter[ , PROTEIN_PADJ := p.adjust(PROTEIN_PVALUE, method = 'BY'), by = PROTEOME_ID ]
    prot_iter[ , SampleSize := i ]
    prot_iter[ , Iteration := j]
    prot_iterations <- rbind(prot_iterations, prot_iter)
  }
}


saveRDS(prot_iterations,  '/data/agprokisch/working/data/proteomics/processed_results/protrider/iterate_MECP2_causal.rds')



