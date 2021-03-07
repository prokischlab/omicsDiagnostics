#'---
#' title: Protrider vs LIMMA      
#' author: smirnovd
#' wb:
#'  input:
#'  - config: 'src/config.R'
#'  - sample_annotation: '`sm config["ANNOTATION"]`'
#'  - protrider_results: '`sm config["PROC_DATA"] + "/protrider/PROTRIDER_results.rds"`'
#'  - limma_results: '`sm config["PROC_DATA"] + "/limma/LIMMA_results.rds"`'
#'  output:
#' output: 
#'   html_document:
#'    code_folding: hide
#'    code_download: TRUE
#'---


# load functions
source(snakemake@input$config)
source("src/functions/LIMMA/limma_functions.R")


# Load sample annotation
sa <- fread(snakemake@input$sample_annotation)
# sa <- fread('/s/project/mitoMultiOmics/multiOMICs_integration/raw_data/proteomics_annotation.tsv')
sa <- sa[USE_FOR_PROTEOMICS_PAPER == T]
sa <- sa[!is.na(KNOWN_MUTATION)]
sa[ , sample_gene := paste0(SAMPLE_ID, "_", KNOWN_MUTATION )]


# Load LIMMA results
limmar <- readRDS(snakemake@input$limma_results) %>% as.data.table()
#limmar <- readRDS("/s/project/mitoMultiOmics/multiOMICs_integration/processed_data/limma/LIMMA_results.rds") %>% as.data.table()
# Remove TASP1 gene. This gene has a bad detection rate and was not detected in this sample. 
# Intensity value was imputed.
limmar <- limmar[!(SAMPLE_ID == "OM30476" & geneID == "TASP1") ]
limmar <- limmar[SAMPLE_ID %in% sa$SAMPLE_ID ]


# Load PROTRIDER results
protr <- readRDS(snakemake@input$protrider_results) %>% as.data.table()
# protr <- readRDS("/s/project/mitoMultiOmics/multiOMICs_integration/processed_data/protrider/PROTRIDER_results.rds") %>% as.data.table()
# Remove TASP1 gene. This gene has a bad detection rate and was not detected in this sample. 
# Intensity value was imputed.
protr <- protr[!(SAMPLE_ID == "OM30476" & geneID == "TASP1") ]
protr <- protr[SAMPLE_ID %in% sa$SAMPLE_ID ]





#combine
prot <- merge(limmar, protr, by = c("SAMPLE_ID", "geneID") , all = T)
prot[, sample_gene := paste0(SAMPLE_ID, "_", geneID )]
prot[, causal_gene := sample_gene %in% sa$sample_gene]

prot[, outlier := "non_outlier"]
prot[PROTEIN_outlier.x == T , outlier := "LIMMA"]
prot[PROTEIN_outlier.y == T , outlier := "PROTRIDER"]
prot[PROTEIN_outlier.x == T  & PROTEIN_outlier.y == T , outlier := "Both"]


outlier_colors <- c(  "non_outlier" = "gray80",
                      "LIMMA" = "#edae49", 
                      "PROTRIDER" = "#FB9A99",  
                      "Both" = "brown2")




range <- c(prot[causal_gene == T ]$PROTEIN_ZSCORE.x, prot[causal_gene == T]$PROTEIN_ZSCORE.y)

#+ fig.width=9, fig.height=9
ggplot(data = prot[causal_gene == T ], aes(PROTEIN_ZSCORE.x, PROTEIN_ZSCORE.y, color= outlier)) +
  geom_point(size = 3) +
  annotate("rect", xmin = -2, xmax = 2, ymin = -2, ymax = 2, alpha = .1, fill = "blue")+
  geom_label_repel(aes(label= geneID , fill = outlier ), color = "gray20", size = 2.8)+ 
  xlab("LIMMA z-score") + 
  ylab("PROTRIDER z-score") +
  geom_abline( color = "grey50") +
  geom_vline(xintercept = 0, color = "grey50") +
  geom_hline(yintercept = 0, color = "grey50") +
  scale_color_manual(breaks = c("non_outlier","LIMMA","PROTRIDER", "Both"), values = outlier_colors)+
  scale_fill_manual(breaks = c("non_outlier","LIMMA","PROTRIDER", "Both"), values = outlier_colors)+
  coord_fixed( xlim = c( min(range , na.rm = T) , max(range , na.rm = T) ),  ylim = c( min(range , na.rm = T),  max(range, na.rm = T)  )) + 
  ggtitle("LIMMA vs PROTRIDER")+
  theme_classic2()+
  theme(plot.title = element_text(size=14, hjust = 0.5, face="bold"),
        axis.title.x= element_text( size=10, margin = NULL,face="bold"),
        axis.title.y= element_text( size=10, margin = NULL,face="bold"),
        axis.text.x= element_text( size=9, margin = NULL,face="bold"),
        axis.text.y= element_text( size=9, margin = NULL,face="bold"))












