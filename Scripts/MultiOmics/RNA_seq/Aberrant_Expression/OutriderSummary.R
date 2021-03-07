#'---
#' title: OUTRIDER Summary
#' author: Michaela Mueller, vyepez
#' wb:
#'  input:
#'   - ods: '`sm config["PROC_DATA"] + "/outrider/ods.Rds"`'
#'   - outrider: '`sm config["PROC_DATA"] + "/outrider/OUTRIDER_results.rds"`'
#'  output:
#' output: 
#'   html_document:
#'    code_folding: hide
#'    code_download: TRUE
#'---

suppressPackageStartupMessages({
    library(OUTRIDER)
    library(vctrs)
    library(SummarizedExperiment)
    library(ggplot2)
    library(data.table)
    library(dplyr)
    library(ggthemes)
    library(readr)
    library(tidyverse)
})


#' ## Read ods object
ods <- readRDS(snakemake@input$ods)

#' ### Aberrant samples
plotAberrantPerSample(ods)


##################
#Load OUTRIDER results 
##################
# rna <- readRDS('/s/project/mitoMultiOmics/multiOMICs_integration/processed_data/outrider/OUTRIDER_results.rds') %>% as.data.table()
rna <- readRDS(snakemake@input$outrider) %>% as.data.table()

#' ### How many samples with at least one outlier gene
rna[ RNA_outlier == T , uniqueN(SAMPLE_ID)]


DT::datatable(rna[RNA_outlier == T], caption = "OUTRIDER results", 
              style = 'bootstrap', filter = 'top', escape = F,
              extensions = c( 'Buttons', 'ColReorder' ), 
              options = list( colReorder = TRUE, dom = 'Bfrtip',
                              buttons = c('copy', 'csv', 'excel', 'pdf', 'print')))







