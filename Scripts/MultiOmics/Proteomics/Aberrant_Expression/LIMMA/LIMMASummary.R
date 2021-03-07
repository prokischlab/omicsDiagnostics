#'---
#' title: LIMMA Summary   
#' author: smirnovd
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
