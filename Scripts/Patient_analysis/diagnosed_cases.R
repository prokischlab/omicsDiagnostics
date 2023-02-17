#'---
#' title: Diagnosed cases 
#' author: Dmitrii Smirnov
#' wb:
#'  input:
#'  - config: 'src/config.R'
#'  - sample_annotation: '`sm config["ANNOTATION"]`'
#'  - patient_omics: '`sm config["PROC_DATA"] + "/integration/patient_omics_full.RDS"`'
#'  output:
#' output: 
#'   html_document:
#'    code_folding: hide
#'    code_download: TRUE
#'---

# Load config
source(snakemake@input$config)


# Load sample annotation
sa <- fread(snakemake@input$sample_annotation)
sa <- sa[USE_FOR_PROTEOMICS_PAPER == T]
sa <- sa[!is.na(KNOWN_MUTATION) & KNOWN_MUTATION != "CONTROL", 1:4]


# Read integrated omics file 
rp <- readRDS(snakemake@input$patient_omics) %>% as.data.table()

# Subset diagnosed cases and candidates
rp <- rp[causal_gene == T]

rpx <- merge(sa, rp, by.x = c("SAMPLE_ID", "KNOWN_MUTATION"), by.y = c("SAMPLE_ID", "geneID"), all = T)

#colnames(rpx)
rpx <- rpx[ , c("SAMPLE_ID", "KNOWN_MUTATION", "gender", "CATEGORY", # "causal_gene",
                                     "gene_class", "outlier_class", "validated", 
                                    "HPO_match", "Semantic_sim", "MDG", "gene_detected",
                                    "RNA_ZSCORE", "RNA_PADJ",  "PROTEIN_ZSCORE", "PROTEIN_PADJ",
                                    "RNA_FC","RNA_PVALUE" , "PROTEIN_FC", "PROTEIN_PVALUE",
                                    "rank_rna","normcounts",
                                    "rank_protein", "PROTEIN_INT",
                                    "WES_avaliable", "HPO_avaliable", 
                                    "RNA_seq_avaliable", "Proteomics_avaliable", "rare", "potential_biallelic")]

rpx <- rpx[ order( CATEGORY)]

#' # Solved cases 
DT::datatable(rpx, 
              caption = "Solved cases", style = 'bootstrap', filter = 'top', escape = F,
              extensions = c( 'Buttons', 'ColReorder' ), 
              options = list( colReorder = TRUE, dom = 'Bfrtip',
                              buttons = c('copy', 'csv', 'excel', 'pdf', 'print')))



