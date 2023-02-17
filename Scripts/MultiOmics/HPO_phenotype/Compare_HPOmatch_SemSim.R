#'---
#' title: Semantic similarity vs HPO match
#' author: smirnovd
#' wb:
#'  input:
#'  - config: 'src/config.R'
#'  - sa: '`sm config["ANNOTATION"]`'
#'  - hpo_match: '`sm config["PROC_DATA"] + "/HPO/Patients_HPO_Gene_mapping.tsv"`'
#'  - semantic_similariy: '`sm config["PROC_DATA"] + "/HPO/Patient_Gene_semantic_similariy.tsv"`'
#'  output:
#'  - phenotype_data: '`sm config["PROC_DATA"] + "/HPO/Patients_phenotype_data.tsv"`'
#' output:
#'   html_document:
#'    code_folding: hide
#'    code_download: TRUE
#'---

source(snakemake@input$config)


# READ ANNOTATION
sa <- fread(snakemake@input$sa)
sa <- sa[USE_FOR_PROTEOMICS_PAPER == T]

# Read HPO match
hpo_match <- fread(snakemake@input$hpo_match)


# Read Semantic similarity scores
SSscores <- fread(snakemake@input$semantic_similariy)


phenotype <- merge(hpo_match, SSscores, by = c("SAMPLE_ID", "geneID"), all = T )
phenotype[ is.na(HPO_match ) , HPO_match := F]
phenotype[ is.na(Semantic_sim ) , Semantic_sim := 0]


write_tsv(phenotype,  snakemake@output$phenotype_data)


ggplot(phenotype, aes(Semantic_sim,  fill = HPO_match))+
  geom_density(alpha= 0.2)+
  theme_classic()+
  geom_vline(xintercept =2, color = "grey70", linetype = "dashed")+
  xlab("Symmetric semantic similarity score") + 
  ylab("Density")+
  scale_fill_ptol()+
  theme(plot.title = element_text(hjust = 0.5),
        axis.title.x= element_text( size=10, margin = NULL,face="bold"),
        axis.title.y= element_text( size=10, margin = NULL,face="bold"),
        axis.text.x = element_text(face="bold",  size=10),
        axis.text.y = element_text(face="bold",  size=10),
        legend.title = element_text(face="bold",  size=10),
        legend.position = c(0.2, 0.8))


ggplot(phenotype, aes(HPO_match ,Semantic_sim,  fill = HPO_match))+
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

