#'---
#' title: PROTRIDER normal Summary
#' author: Dmitrii Smirnov
#' wb:
#'  input:
#'   - protrider_resultsu: '`sm config["PROC_DATA"] + "/protrider/PROTRIDER_results_gaus.rds"`'
#'  output:
#' output: 
#'   html_document:
#'    code_folding: hide
#'    code_download: TRUE
#'---

source("src/functions/LIMMA/limma_functions.R")
library(data.table)
library(DT)

prot <- readRDS(snakemake@input$protrider_resultsu) %>% as.data.table()

 
# Total number of outliers:
uniqueN(prot[PROTEIN_outlier == TRUE])
paste("Underexpression:", uniqueN(prot[PROTEIN_outlier == TRUE & PROTEIN_ZSCORE < 0]),
      ", Overexpression:", uniqueN(prot[PROTEIN_outlier == TRUE & PROTEIN_ZSCORE > 0]))

# Aberrant proteins per sample:
os <- prot[PROTEIN_outlier == TRUE, .N, by = SAMPLE_ID]
paste("Median:", median(os$N))
plotAberrantProteinPerSample(os)

# How many samples with at least one aberrant protein:
prot[PROTEIN_outlier == TRUE, uniqueN(SAMPLE_ID)]


DT::datatable(prot[PROTEIN_outlier == TRUE], caption = "PROTRIDER results", 
              style = 'bootstrap', filter = 'top', escape = FALSE,
              extensions = c('Buttons', 'ColReorder'),
              options = list(colReorder = TRUE, dom = 'Bfrtip',
                             buttons = c('copy', 'csv', 'excel', 'pdf', 'print')))
