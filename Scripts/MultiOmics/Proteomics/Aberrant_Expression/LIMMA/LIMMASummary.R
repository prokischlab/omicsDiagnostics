#'---
#' title: LIMMA Summary   
#' author: Dmitrii Smirnov
#' wb:
#'  input:
#'  - results: '`sm config["PROC_DATA"] + "/limma/LIMMA_results.rds"`'
#'  output:
#' output: 
#'   html_document:
#'    code_folding: hide
#'    code_download: TRUE
#'---

# load functions
source("src/config.R")
source("src/functions/LIMMA/limma_functions.R")


# Results
prot <- readRDS(snakemake@input$results) %>% as.data.table()


# Total number of outliers:
uniqueN(prot[PROTEIN_outlier == TRUE])
paste("Underexpression:", uniqueN(prot[PROTEIN_outlier == TRUE & PROTEIN_ZSCORE < 0]),
      ", Overexpression:", uniqueN(prot[PROTEIN_outlier == TRUE & PROTEIN_ZSCORE > 0]))

# Aberrant proteins per sample:
os <- prot[PROTEIN_outlier == TRUE, .N, by = SAMPLE_ID]
paste("Median:", median(os$N))
plotAberrantProteinPerSample(os)


#' ### Aberrant proteins per sample
plotAberrantProteinPerSample(prot[ PROTEIN_outlier == T , .N, by = c('SAMPLE_ID')])


#' ### How many samples with at least one outlier gene
prot[PROTEIN_outlier == T, uniqueN(SAMPLE_ID)]


#+echo=F
DT::datatable(prot[PROTEIN_outlier == T], caption = "LIMMA results", 
              style = 'bootstrap', filter = 'top', escape = F,
              extensions = c( 'Buttons', 'ColReorder' ), 
              options = list( colReorder = TRUE, dom = 'Bfrtip',
                              buttons = c('copy', 'csv', 'excel', 'pdf', 'print')))
