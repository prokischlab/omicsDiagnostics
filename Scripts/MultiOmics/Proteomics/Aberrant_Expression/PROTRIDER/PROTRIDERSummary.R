#'---
#' title: PROTRIDER Summary
#' author: Dmitrii Smirnov
#' wb:
#'  input:
#'  - protrider_results: '`sm config["PROC_DATA"] + "/protrider/PROTRIDER_results.rds"`'
#'  output:
#' output: 
#'   html_document:
#'    code_folding: hide
#'    code_download: TRUE
#'---


# source("src/config.R")
source("src/functions/LIMMA/limma_functions.R")


#' ## Results
prot <- readRDS(snakemake@input$protrider_results) %>% as.data.table()


#' Total number of outliers: 
uniqueN(prot[ PROTEIN_outlier == T])
paste("Underexpression:",  uniqueN(prot[ PROTEIN_outlier == T  & PROTEIN_ZSCORE <0 ]), ", Overexpression:", uniqueN(prot[ PROTEIN_outlier == T & PROTEIN_ZSCORE > 0 ])  )

#' ### Aberrant proteins per sample
#' 
os <- prot[ PROTEIN_outlier == T , .N, by = SAMPLE_ID ]
paste("Median:",  median(os$N) )
plotAberrantProteinPerSample(os)



#' ### How many samples with at least one gene
prot[PROTEIN_outlier == T, uniqueN(SAMPLE_ID)]


#+echo=F
DT::datatable(prot[PROTEIN_outlier == T], caption = "PROTRIDER results", 
              style = 'bootstrap', filter = 'top', escape = F,
              extensions = c( 'Buttons', 'ColReorder' ), 
              options = list( colReorder = TRUE, dom = 'Bfrtip',
                              buttons = c('copy', 'csv', 'excel', 'pdf', 'print')))









