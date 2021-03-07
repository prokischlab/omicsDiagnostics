#'---
#' title: Figure 2a RNA vs protein integration all individuals
#' author: smirnovd
#' wb:
#'  input:
#'  - config: 'src/config.R'
#'  - patient_omics: '`sm config["PROC_DATA"] + "/integration/patient_omics.RDS"`'
#'  output:
#'  - fig: '`sm config["FIGURE_DIR"] + "/Fig2_d.pdf"`'
#' output: 
#'   html_document:
#'    code_folding: hide
#'    code_download: TRUE
#'---

source(snakemake@input$config)


# Read integrated omics file 
# rp <- readRDS("/s/project/mitoMultiOmics/multiOMICs_integration/processed_data/integration/patient_omics.RDS") %>% as.data.table()
rp <- readRDS(snakemake@input$patient_omics) %>% as.data.table()


rp <- rp[gene_class != "no data"]

#+ fig.width=8, fig.height=8
#range <- c(rp[outlier_class != 'non_outlier' ]$RNA_ZSCORE, rp[outlier_class != 'non_outlier' ]$PROTEIN_ZSCORE)

fig <- ggplot(data = rp[gene_detected == "RNA and protein detected"], aes(RNA_ZSCORE, PROTEIN_ZSCORE, color= outlier_class, shape= gene_class)) +
  geom_hex(data = rp[gene_detected == "RNA and protein detected" & outlier_class == 'non_outlier' ], 
           aes(RNA_ZSCORE, PROTEIN_ZSCORE, fill= outlier_class, color= outlier_class), bins=100, show.legend = F) +
  
  geom_point(data = rp[gene_detected == "RNA and protein detected" & outlier_class != 'non_outlier' & gene_class == "no rare variant"], 
             aes(RNA_ZSCORE, PROTEIN_ZSCORE, color = outlier_class, 
                 shape= gene_class,  size = Semantic_sim ), show.legend = F) +
  
  geom_point(data = rp[gene_detected == "RNA and protein detected" & outlier_class != 'non_outlier' & gene_class != "no rare variant" ], 
             aes(RNA_ZSCORE, PROTEIN_ZSCORE, fill= outlier_class, 
                 shape= gene_class, size = Semantic_sim ), color= "grey20") +
  xlab("RNA zScore") + 
  ylab("Protein zScore") +
  geom_vline(xintercept = 0, color = "grey50") +
  geom_hline(yintercept = 0, color = "grey50") +
  scale_color_manual(breaks = c("non_outlier","RNA_outlier","Protein_outlier", "RNA_Protein_outlier"),
                     values = outlier_colors)+
  scale_fill_manual(breaks = c("non_outlier","RNA_outlier","Protein_outlier", "RNA_Protein_outlier"),
                    values = outlier_colors)+
  
  scale_shape_manual(values = c("no rare variant" = 20,  
                                "1 rare variant" = 24, 
                                "rare pot. biallelic variants" = 22) ) +
  
  coord_fixed( xlim = c( -11 ,11),  ylim = c( -11,  11  )) + 
  theme_bw()+
  scale_size(range = c(2.5,6))+
  labs(color='Outlier class',
       fill='Outlier class',
       size='Semantic similarity',
       shape='Variant information') +
  guides(color=guide_legend(nrow=4,byrow=TRUE),
         fill=guide_legend(nrow=4,byrow=TRUE),
         size=guide_legend(nrow=6,byrow=TRUE),
         shape=guide_legend(nrow=3,byrow=TRUE))+
  theme(axis.title.x= element_text( size=10, margin = NULL,face="bold"),
        axis.title.y= element_text( size=10, margin = NULL,face="bold"),
        axis.text.x = element_text(face="bold",  size=10),
        axis.text.y = element_text(face="bold",  size=10),
        legend.position="bottom",
        plot.margin = margin(0, 0, 0, 0, "cm"))
fig

pdf(snakemake@output$fig, # "/s/project/mitoMultiOmics/multiOMICs_integration/Figures/Fig2_a.pdf",  
    width = 7, height = 7,  useDingbats=FALSE )
print(fig) 
dev.off()

