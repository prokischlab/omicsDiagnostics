#'---
#' title: Figure 2c MORC2 dominant case example
#' author: smirnovd
#' wb:
#'  input:
#'  - config: 'src/config.R'
#'  - patient_omics: '`sm config["PROC_DATA"] + "/integration/patient_omics.RDS"`'
#'  - sample_annotation: '`sm config["ANNOTATION"]`'
#'  - patient_organs: '`sm config["PROC_DATA"] + "/HPO/Patients_affected_organs.tsv"`'
#'  output:
#'  - fig_human: '`sm config["FIGURE_DIR"] + "/Fig2_c_human.pdf"`'
#'  - fig: '`sm config["FIGURE_DIR"] + "/Fig2_c.pdf"`'
#' output: 
#'   html_document:
#'    code_folding: hide
#'    code_download: TRUE
#'---

dev.off()
source(snakemake@input$config)


# Add necessaries for gene structures
suppressPackageStartupMessages({
  stringsAsFactors=F
  library(Gviz)
  library(GenomicRanges)
  library(ggbio)
  library(EnsDb.Hsapiens.v75)
})



## Making a short cut to ensembl 
edb <- EnsDb.Hsapiens.v75


# Load sample annotation
# sa <- fread('/s/project/mitoMultiOmics/multiOMICs_integration/raw_data/proteomics_annotation.tsv')
sa <- fread(snakemake@input$sample_annotation)
sa <- sa[USE_FOR_PROTEOMICS_PAPER == T]

# Read integrated omics file 
# rp <- readRDS("/s/project/mitoMultiOmics/multiOMICs_integration/processed_data/integration/patient_omics.RDS") %>% as.data.table()
rp <- readRDS(snakemake@input$patient_omics) %>% as.data.table()


# Load patient's affected organs
patients <- fread(snakemake@input$patient_organs)
# patients <- fread('/s/project/mitoMultiOmics/multiOMICs_integration/processed_data/HPO/Patients_affected_organs.tsv')


######################################################
#+ fig.width=3, fig.height=5

fig_human <- gganatogram(data=  patients[SAMPLE_ID == "OM27390"], 
                           fillOutline='white', organism='human', 
                           sex= unique(patients[SAMPLE_ID == "OM27390"]$sex), 
                           fill="colour")+ 
  theme_void() + 
  theme(plot.title = element_text(margin = NULL,face="bold", size = 9, hjust = 0.5),
        plot.margin = margin(0, 0, 0, 0, "cm"))
fig_human

pdf(snakemake@output$fig_human, # "/s/project/mitoMultiOmics/multiOMICs_integration/Figures/Fig2_c_human.pdf",  
    width = 3, height =5,  useDingbats=FALSE )
print(fig_human) 
dev.off()


samp <- rp[ SAMPLE_ID == "OM27390" &  (outlier_class != "non_outlier" | causal_gene == T ) ]
fig <- ggplot(samp[gene_detected == "RNA and protein detected" ], aes(RNA_ZSCORE, PROTEIN_ZSCORE)) +
  geom_vline(xintercept = 0, color = "grey50") +
  geom_hline(yintercept = 0, color = "grey50") +
  
  geom_point(aes(color= outlier_class, shape= gene_class, size = Semantic_sim )) + 
  
  xlab("RNA zScore") + 
  ylab("Protein zScore") +
  scale_color_manual(breaks = c("non_outlier","RNA_outlier","Protein_outlier", "RNA_Protein_outlier"),
                     values = outlier_colors) +
  
  scale_shape_manual(values = gene_shapes) +
  
  scale_size(range = c(2.5, 6)) +
  
  geom_text_repel(data= samp[gene_detected == "RNA and protein detected"  & gene_class == "no rare variant" ], 
                  mapping=aes(RNA_ZSCORE, PROTEIN_ZSCORE, label= geneID), colour= "grey20", #casual_gene
                  size= 2.5, show.legend = F)+ 
  
  geom_text_repel(data= samp[gene_detected == "RNA and protein detected"  & gene_class == "1 rare variant" & ( causal_gene ==F | is.na(causal_gene)) ], 
                  mapping=aes(RNA_ZSCORE, PROTEIN_ZSCORE, label= geneID), colour= "black", 
                  size= 3.1, show.legend = F) + 
  
  geom_label_repel(data= samp[gene_detected == "RNA and protein detected" & gene_class == "rare pot. biallelic variants" | causal_gene == T], 
                   mapping=aes(RNA_ZSCORE, PROTEIN_ZSCORE, label= geneID ), 
                   box.padding = unit(0.3, "lines"), colour= "black", 
                   size= 4, show.legend = F ) + 
  
  coord_fixed( xlim = c( -8, 8 ),  ylim = c( -8,  8 )) +  
  theme_bw() +
  
  theme( axis.title.x= element_text( size= 10, margin = NULL,face="bold"),
         axis.title.y= element_text( size= 10, margin = NULL,face="bold"),
         legend.title = element_blank(),
         axis.text.x = element_text(face="bold",  size=10),
         axis.text.y = element_text(face="bold",  size=10),
         legend.position = "none" )  # c(0.2, 0.8)

#+ fig.width=4, fig.height=4
fig

pdf(snakemake@output$fig, # "/s/project/mitoMultiOmics/multiOMICs_integration/Figures/Fig2_c.pdf",  
    width = 4, height =4,  useDingbats=FALSE )
print(fig) 
dev.off()



#' Gene structure and variant 
#+ fig.width=5, fig.height=2

autoplot(edb, ~ gene_name == "MORC2", stat = "reduce" ,  names.expr="MORC2") + 
  theme_void() +  scale_x_reverse() + 
  geom_vline( xintercept =  31345792 , colour = "darkblue")
dev.off()


unloadNamespace("ggbio")



