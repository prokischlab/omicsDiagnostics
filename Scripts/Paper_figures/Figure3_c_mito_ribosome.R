#'---
#' title: Figure 3c MRPL38 mitochondrial ribosome  
#' author: Dmitrii Smirnov
#' wb:
#'  input:
#'  - config: 'src/config.R'
#'  - patient_omics: '`sm config["PROC_DATA"] + "/integration/patient_omics_full.RDS"`'
#'  - corum: '`sm config["DATASETS"] + "/CORUM.tsv"`'
#'  output:
#'  - fig: '`sm config["FIGURE_DIR"] + "/Fig3_c.pdf"`'
#' output: 
#'   html_document:
#'    code_folding: hide
#'    code_download: TRUE
#'---

source(snakemake@input$config)

# Read integrated omics file 
rp <- readRDS(snakemake@input$patient_omics) %>% as.data.table()




# Read CORUM complexes
corum <- fread(snakemake@input$corum)

## Get Corum MRPL/S complexes
corum <- corum[ COMPLEX %in% c("28S ribosomal subunit, mitochondrial", "39S ribosomal subunit, mitochondrial")]

case_MRPL38 <- rp[SAMPLE_ID == "OM57837"] 
case_MRPL38 <- case_MRPL38[!is.na(PROTEIN_LOG2FC )]

case_MRPL38[, small_sub := geneID %in% corum[COMPLEX == "28S ribosomal subunit, mitochondrial"]$geneID]
case_MRPL38[, large_sub := geneID %in% corum[COMPLEX == "39S ribosomal subunit, mitochondrial"]$geneID]


# Small subunit highlighted in volcano plot
ss <- ggplot(case_MRPL38, aes(PROTEIN_LOG2FC, -log10(PROTEIN_PVALUE) ))+
  geom_vline(xintercept = 0, color = "grey30") +
  geom_point( color = "gray80") +
  geom_point( data = case_MRPL38[small_sub == T ],
              aes(PROTEIN_LOG2FC, -log10(PROTEIN_PVALUE)), color = "#edae49" ) + 
  xlab("log2 fold change") + 
  ylab("-log10(p-value)") +
  ggtitle("small ribosomal proteins")+
  scale_x_continuous(limits = c( -max( abs(case_MRPL38$PROTEIN_LOG2FC) , na.rm = T), 
                                 max( abs(case_MRPL38$PROTEIN_LOG2FC) , na.rm = T) ))+  
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5, size=12, color = "#edae49", face="bold"),
        axis.title.x= element_text( size=10, margin = NULL,face="bold"),
        axis.title.y= element_text( size=10, margin = NULL,face="bold"),
        axis.text.x = element_text(face="bold",  size=10),
        axis.text.y = element_text(face="bold",  size=10),
        legend.title = element_blank() )

# Large subunit highlighted in volcano plot
ls <- ggplot(case_MRPL38, aes(PROTEIN_LOG2FC, -log10(PROTEIN_PVALUE) ))+
  geom_vline(xintercept = 0, color = "grey30") +
  geom_point(color = "gray80") +
  geom_point( data = case_MRPL38[large_sub == T ],
              aes(PROTEIN_LOG2FC, -log10(PROTEIN_PVALUE)), color = "#00798c" ) + 
  geom_point(data = case_MRPL38[ geneID == "MRPL38"], 
             aes(PROTEIN_LOG2FC, -log10(PROTEIN_PVALUE) ), color = "#FB9A99" ) +
  geom_text_repel(data = case_MRPL38[ geneID == "MRPL38"], aes(PROTEIN_LOG2FC, -log10(PROTEIN_PVALUE), label =  geneID), size = 3,fontface ="bold" ) +
  xlab("log2 fold change") + 
  ylab("-log10(p-value)") +
  ggtitle("large ribosomal proteins")+
  
  scale_x_continuous(limits = c( -max( abs(case_MRPL38$PROTEIN_LOG2FC) , na.rm = T), 
                                 max( abs(case_MRPL38$PROTEIN_LOG2FC) , na.rm = T) ))+ 
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5, size=12, colour = "#00798c", face="bold"),
        axis.title.x= element_text( size=10, margin = NULL,face="bold"),
        axis.title.y= element_text( size=10, margin = NULL,face="bold"),
        axis.text.x = element_text(face="bold",  size=10),
        axis.text.y = element_text(face="bold",  size=10),
        legend.title = element_blank() )

#+ fig.width=5, fig.height=7
fig <- ss / ls
fig


pdf(snakemake@output$fig, 
    width = 5, height =7,  useDingbats=FALSE )
print(fig) 
dev.off()

