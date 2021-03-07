#'---
#' title: Figure 1d VUS validation
#' author: smirnovd
#' wb:
#'  input:
#'  - config: 'src/config.R'
#'  - patient_omics: '`sm config["PROC_DATA"] + "/integration/patient_omics.RDS"`'
#'  - sample_annotation: '`sm config["ANNOTATION"]`'
#'  output:
#'  - fig: '`sm config["FIGURE_DIR"] + "/Fig1_d.pdf"`'
#' output: 
#'   html_document:
#'    code_folding: hide
#'    code_download: TRUE
#'---


# load functions
source(snakemake@input$config)

dev.off()


##########################################################
# Read integrated omics file 
# rp <- readRDS("/s/project/mitoMultiOmics/multiOMICs_integration/processed_data/integration/patient_omics.RDS") %>% as.data.table()
rp <- readRDS(snakemake@input$patient_omics) %>% as.data.table()


# Load sample annotation
# sa <- fread('/s/project/mitoMultiOmics/multiOMICs_integration/raw_data/proteomics_annotation.tsv')
sa <- fread(snakemake@input$sample_annotation)
sa <- sa[USE_FOR_PROTEOMICS_PAPER == T]
sa$sample_gene <- paste0(sa$SAMPLE_ID, "_", sa$KNOWN_MUTATION)
sa <- sa[order(CATEGORY, KNOWN_MUTATION)] 
sax <- sa


# Subset cases with validated VUS (IIa)
sa <- sa[CATEGORY %in% c("IIa") ]
sa$sample_gene <- paste0(sa$SAMPLE_ID, "_", sa$KNOWN_MUTATION)


# Subset data for causal genes in all samples
rp <- rp[ geneID %in% sa$KNOWN_MUTATION]
rp <- rp[!duplicated(rp)]
rp$sample_gene <- paste0(rp$SAMPLE_ID, "_", rp$geneID)


# Re-annotate causal gene
rp$causal_gene <- NULL
rp[, causal_gene :=  sample_gene %in% sa$sample_gene]

# Some genes reappear in discovery cohort 
rp[, another_cohort := sample_gene %in% sax[CATEGORY %in% c("I",  "III") ]$sample_gene ]


rp[causal_gene == T & is.na(PROTEIN_ZSCORE), PROTEIN_ZSCORE := -Inf]
rp[causal_gene == T & is.na(RNA_ZSCORE), RNA_ZSCORE := -Inf]


# Sort by causal gene Z score
rp_c <- rp[causal_gene == T ]
rp$geneID <- factor(rp$geneID, unique(rp_c[ order(PROTEIN_ZSCORE)]$geneID))

#####################################################################


#+ fig.width=4, fig.height=3
fig <- ggplot(rp, aes(geneID,  PROTEIN_ZSCORE ))+
  geom_hline(yintercept = 0, color = "grey50") +
  geom_quasirandom(data = rp[causal_gene != T ] , aes(geneID,  PROTEIN_ZSCORE),color = "grey80" , size=0.5) +
  stat_summary(fun.data = mean_sdl, fun.args = list(mult=1), geom="errorbar", color="black", width=0.2) +
  stat_summary(fun = mean, geom="point", color="black") +
  ylab("Protein zScore") + 
  scale_y_continuous(limits = c( -10,  5 ))+ 
  geom_point( data = rp[causal_gene == T ], aes(geneID,  PROTEIN_ZSCORE), colour = "#FB9A99" , size = 2) +
  geom_point( data = rp[another_cohort == T ], aes(geneID,  PROTEIN_ZSCORE ), size = 1, colour = "darkgreen", shape = 8)+
  theme_classic()+
  theme( legend.title = element_blank(),
         axis.title.x = element_blank(),
         plot.title = element_text(face="bold",  size=14, hjust = 0.5),
         axis.title.y = element_text(face="bold",  size=10),
         axis.text.x = element_text(face="bold",  size=9, angle = 45, hjust = 1),
         axis.text.y = element_text(face="bold",  size=9),
         legend.position = "none")
fig


pdf(snakemake@output$fig, # "/s/project/mitoMultiOmics/multiOMICs_integration/Figures/Fig1_d.pdf",  
    width = 4, height = 3,  useDingbats=FALSE )
print(fig) 
dev.off()









