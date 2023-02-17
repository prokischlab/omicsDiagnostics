#'---
#' title: Supplementary Fig 2d Semantic similarity vs HPO_match
#' author: Dmitrii Smirnov
#' wb:
#'  input: 
#'  - config: 'src/config.R'
#'  - var_hpo: '`sm config["RAW_DATA"] + "/patient_variant_hpo_data.tsv"`'
#'  output:
#'  - fig: '`sm config["FIGURE_DIR"] + "/Supplementary_figures/S_Fig2_d.pdf"`'
#' output: 
#'   html_document:
#'    code_folding: hide
#'    code_download: TRUE
#'---


# load functions
source(snakemake@input$config)


# Load phenotype information for all cases
var_hpo <- fread(snakemake@input$var_hpo)
var_hpo <- var_hpo[!is.na(HPO_match)]

s_fig <- ggplot(var_hpo, aes(HPO_match ,Semantic_sim,  fill = HPO_match))+
  geom_boxplot(alpha= 0.7)+ 
  theme_classic()+
  geom_hline(yintercept =2, color = "grey70", linetype = "dashed")+
  scale_fill_ptol()+
  ylab("Symmetric semantic similarity score") + 
  xlab("HPO match")+
  theme(plot.title = element_text(hjust = 0.5),
        axis.title.x= element_text( size=10, margin = NULL,face="bold"),
        axis.title.y= element_text( size=10, margin = NULL,face="bold"),
        axis.text.x = element_text(face="bold",  size=10),
        axis.text.y = element_text(face="bold",  size=10),
        legend.title = element_blank(),
        legend.position = "none")

#+ fig.width=6, fig.height=6
s_fig


pdf(snakemake@output$fig, 
    width = 6, height =6,  useDingbats=FALSE )
print(s_fig) 
dev.off()





