#'---
#' title: Figure 1b rna and proteome coverage
#' author: smirnovd
#' wb:
#'  input: 
#'  - protein_coverage: '`sm config["PROC_DATA"] + "/integration/protein_coverage.tsv"`'
#'  - rna_coverage: '`sm config["PROC_DATA"] + "/integration/rna_coverage.tsv"`'
#'  output:
#'  - fig: '`sm config["FIGURE_DIR"] + "/Fig1_b.pdf"`'
#' output: 
#'   html_document:
#'    code_folding: hide
#'    code_download: TRUE
#'---

# Load config
source('src/config.R')


# rna_c <- fread('/s/project/mitoMultiOmics/multiOMICs_integration/processed_data/integration/rna_coverage.tsv')
rna_c <- fread(snakemake@input$rna_coverage)
rna_c[, type:= "RNAseq"]


# prot_c <- fread('/s/project/mitoMultiOmics/multiOMICs_integration/processed_data/integration/protein_coverage.tsv')
prot_c <- fread(snakemake@input$protein_coverage)
prot_c[, type:= "Proteomics"]




omics <- rbind(rna_c, prot_c)
omics[ , DETECTED_omic := paste0(DETECTED, "_", type)]

omics$type <- factor(omics$type , levels = c("RNAseq", "Proteomics"))  

              
omics$DETECTED_omic <- factor(omics$DETECTED_omic , 
                                 levels = c( "ONCE_RNAseq",
                                             "HALF_RNAseq", 
                                             "ALL_RNAseq",
                                             
                                             "ONCE_Proteomics",
                                             "HALF_Proteomics", 
                                             "ALL_Proteomics" 
                                            )) 
omics$dis_n <- factor(omics$dis_n, levels = c("Protein coding\n(20336)", "MITO\n(388)", "Neuromuscular\n(132)", "Neurology\n(284)", "OMIM\n(4270)" ))


omics[DISEASE ==  "Protein coding", protein_group := paste0("protein coding\ngenes\n", "(n = ", total , ")")]
omics[DISEASE ==  "MITO", protein_group := paste0("Mito disease\ngenes\n", "(n = ", total , ")")]
omics[DISEASE ==  "Neuromuscular", protein_group := paste0("Neuromuscular\ndisease genes\n", "(n = ", total , ")")]
omics[DISEASE ==  "Neurology", protein_group := paste0("Neurology\ndisease genes\n", "(n = ", total , ")")]
omics[DISEASE ==  "OMIM", protein_group := paste0("OMIM\ngenes\n", "(n = ", total , ")")]
unique(omics$protein_group)

omics$protein_group <- factor(omics$protein_group, 
                              levels = c("protein coding\ngenes\n(n = 20336)", 
                                         "Mito disease\ngenes\n(n = 388)", 
                                         "Neuromuscular\ndisease genes\n(n = 132)", 
                                         "Neurology\ndisease genes\n(n = 284)", 
                                         "OMIM\ngenes\n(n = 4270)" ))


coverage_colors <- c( "ONCE_RNAseq"= "lightblue",
                      "HALF_RNAseq"= "#96C7E8", 
                      "ALL_RNAseq" = "#2B6EB3",
                      
                      "ONCE_Proteomics"= "pink",
                      "HALF_Proteomics"= "#F06B6E", 
                      "ALL_Proteomics" = "#DB262E")

# write_tsv(omics, "RNAseq_proteomics_coverage.tsv")

#' # RNA-seq and proteomics coverage in fibroblasts
#+ fig.width=7, fig.height=4

fig <- ggplot(omics, aes(type, prop)) + 
  geom_col( aes(fill = DETECTED_omic)) +  #stat= 'identity',
  geom_hline(yintercept = 0, colour = "black") + 
  scale_y_continuous( labels=scales::percent) +
  scale_fill_manual(values = coverage_colors  ) + 
  labs( y = "fraction covered")+
  theme_classic()+
  facet_wrap(~protein_group,  ncol= 5, strip.position = "bottom")+
  theme(  
        axis.title.x = element_blank() ,
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        
        axis.title.y = element_text(face="bold", size=12) ,
        axis.text.y = element_text(size=10, face="bold") ,
        
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        
        legend.position="none",   
        legend.title = element_blank(),
        
        strip.background = element_blank(),
        strip.text = element_text(face="bold", size=11),
        plot.background = element_rect( fill = "white")) 


fig



pdf(snakemake@output$fig, # "/s/project/mitoMultiOmics/multiOMICs_integration/Figures/Fig1_b.pdf",  
    width = 7, height =4,  useDingbats=FALSE )
print(fig) 
dev.off()



#' # Values
DT::datatable(omics[ , c("type", "DISEASE", "total",  "N", "prop")], 
              caption = "", style = 'bootstrap', filter = 'top', escape = F,
              extensions = c( 'Buttons', 'ColReorder' ), 
              options = list( colReorder = TRUE, dom = 'Bfrtip',
                              buttons = c('copy', 'csv', 'excel', 'pdf', 'print')))


