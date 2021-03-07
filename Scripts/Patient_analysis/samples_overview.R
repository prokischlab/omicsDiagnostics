#'---
#' title: Samples overview  
#' author: smirnovd
#' wb:
#'  input:
#'  - sample_annotation: '`sm config["ANNOTATION"]`'
#'  output:
#' output: 
#'   html_document:
#'    code_folding: hide
#'    code_download: TRUE
#'---


# Load config
source("src/config.R")




# Load sample annotation
# sa <- fread('/s/project/mitoMultiOmics/multiOMICs_integration/raw_data/proteomics_annotation.tsv')
sa <- fread(snakemake@input$sample_annotation)



#' # Proteomics sample annotation
DT::datatable(sa, 
              caption = "Proteomics annotation", style = 'bootstrap', filter = 'top', escape = F,
              extensions = c( 'Buttons', 'ColReorder' ), 
              options = list( colReorder = TRUE, dom = 'Bfrtip',
                              buttons = c('copy', 'csv', 'excel', 'pdf', 'print')))



sa <- sa[USE_FOR_PROTEOMICS_PAPER == T]

#' ## Cases per category
print(paste("Total:", nrow(sa)))
print(paste("Healthy:", nrow(sa[AFFECTED == F])))
print(paste("Solved:", nrow(sa[ CATEGORY %in%  c("I", "IIa", "III" )])))


sa[AFFECTED == F , CATEGORY:= "Healthy"]
sa[CATEGORY == "I" , CATEGORY:= "Category I- published variants / Clinvar pathogenic"]
sa[CATEGORY == "I.m" , CATEGORY:= "Category Im - mtDNA cases"]
cat2 <- nrow(sa[CATEGORY %in% c("IIa", "IIb", "IIc") ])
sa[CATEGORY == "IIa" , CATEGORY:= "Category IIa - validated VUS"]
sa[CATEGORY == "IIb" , CATEGORY:= "Category IIb - not validated VUS"]
sa[CATEGORY == "IIc" , CATEGORY:= "Category IIc - VUS not detected by proteomics"]
sa[CATEGORY == "III" , CATEGORY:= "Category III - new discoveries"]


os <- sa[!is.na(CATEGORY) , .N , by = CATEGORY]
os <- os[order(CATEGORY)]
print(paste("Category II - VUS:", cat2 ))

#'
DT::datatable(os, 
              caption = "Cases per category", style = 'bootstrap', filter = 'top', escape = F,
              extensions = c( 'Buttons', 'ColReorder' ), 
              options = list( colReorder = TRUE, dom = 'Bfrtip',
                              buttons = c('copy', 'csv', 'excel', 'pdf', 'print')))


