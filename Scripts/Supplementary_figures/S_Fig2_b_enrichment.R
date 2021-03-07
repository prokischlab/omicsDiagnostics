#'---
#' title: Supplementary Fig 2b Enrichment of rare variants in outliers
#' author: smirnovd
#' wb:
#'  input: 
#'  - config: 'src/config.R'
#'  - enrichment_rare_pb: '`sm config["PROC_DATA"] + "/variant_tables/enrichment_rare_biallelic.tsv"`'
#'  output:
#'  - fig: '`sm config["FIGURE_DIR"] + "/Supplementary_figures/S_Fig2_b.pdf"`'
#' output: 
#'   html_document:
#'    code_folding: hide
#'    code_download: TRUE
#'---


# load config
source(snakemake@input$config)

# Read enrichment results
# enrichments <- fread("/s/project/mitoMultiOmics/multiOMICs_integration/processed_data/variant_tables/enrichment_rare_biallelic.tsv")
enrichments <- fread(snakemake@input$enrichment_rare_pb)
nodelist <- levels(enrichments$Cat)
enrichments$Cat <- factor(enrichments$Cat, levels= c("No rare variant", "Rare variant", "Potential biallelic\n rare variant"))





#+ fig.width=6, fig.height=6
s_fig <- ggplot(data = enrichments, aes(x = Cat, y = Estim)) +
  geom_pointrange(aes(x = Cat, ymin = ci_left, ymax = ci_right, colour = outlier_class), position = position_dodge(width = 1)) +
  geom_hline(yintercept = 0) +
  geom_vline(color = "gray80", linetype = "dashed", xintercept = which(nodelist=="No rare variant")-0.5) +
  geom_vline(color = "gray80", linetype = "dashed", xintercept = which(nodelist=="Rare variant")-0.5) +
  theme_bw() + 
  xlab('') + 
  ylab('Log odds ratio')  + 
  scale_alpha(guide = 'none')+
  scale_color_manual(values = outlier_colors)+ 
  coord_flip(ylim = c(min(enrichments[enrichments$Pval <0.05, 'ci_left' ] , na.rm = T),  
                      max(enrichments[enrichments$Pval <0.05, "ci_right" ], na.rm = T)))+
  theme(plot.title = element_text(hjust = 0.5),
        axis.title.x= element_text( size=10, margin = NULL,face="bold"),
        axis.title.y= element_text( size=10, margin = NULL,face="bold"),
        axis.text.x = element_text(face="bold",  size=10),
        axis.text.y = element_text(face="bold",  size=10),
        legend.title = element_blank())


#+ fig.width=6, fig.height=6
s_fig


pdf(snakemake@output$fig, # "/s/project/mitoMultiOmics/multiOMICs_integration/Figures/Supplementary_figures/S_Fig2_b.pdf",  
    width = 6, height =6,  useDingbats=FALSE )
print(s_fig) 
dev.off()





