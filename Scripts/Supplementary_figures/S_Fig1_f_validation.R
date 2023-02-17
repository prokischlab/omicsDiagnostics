#'---
#' title: Supplementary Fig 1f validation of positive controls and VUS
#' author: Dmitrii Smirnov
#' wb:
#'  input: 
#'  - config: 'src/config.R'
#'  - patient_omics: '`sm config["PROC_DATA"] + "/integration/patient_omics.RDS"`'
#'  - sample_annotation: '`sm config["ANNOTATION"]`'
#'  output:
#'  - fig: '`sm config["FIGURE_DIR"] + "/Supplementary_figures/S_Fig1_f.pdf"`'
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
sa <- fread(snakemake@input$sample_annotation)
sa <- sa[USE_FOR_PROTEOMICS_PAPER == T]
sa$sample_gene <- paste0(sa$SAMPLE_ID, "_", sa$KNOWN_MUTATION)
sa <- sa[order(CATEGORY, KNOWN_MUTATION)] 


# Remove UQCRFS1 as soon as it was not detected by proteomics 
# and RNA-seq is not avaliable for the case with known pathogenic vaiant 
sa <- sa[KNOWN_MUTATION != "UQCRFS1"]


# Subset cases with pathogenic variants (I) and VUS (II)
sa <- sa[CATEGORY %in% c( "I", "IIa", "IIb", "IIc") ]
sa[CATEGORY %in% c("IIa", "IIb", "IIc"), CATEGORY := "II"]
sa$CATEGORY <- factor(sa$CATEGORY, c( "I", "II"))


# Subset data for causal genes in all samples
rp <- rp[ geneID %in% sa$KNOWN_MUTATION]



# ACAD9 appears in both categories, duplicate the entry
sa[KNOWN_MUTATION == "ACAD9" & SAMPLE_ID == "OM65708", KNOWN_MUTATION := " ACAD9"]
sa$sample_gene <- paste0(sa$SAMPLE_ID, "_", sa$KNOWN_MUTATION)

rp_acad9 <- rp[geneID == "ACAD9" & SAMPLE_ID != "OM65728"]
rp_acad9$geneID <- paste0(" ", rp_acad9$geneID)
rp <- rp[ !(geneID == "ACAD9" & SAMPLE_ID == "OM65708") ]
rp <- rbind(rp, rp_acad9 )
rp$sample_gene <- paste0(rp$SAMPLE_ID, "_", rp$geneID)
rm(rp_acad9)


# Annotate with category
gene_cat <- sa[ , c("KNOWN_MUTATION", "CATEGORY")]
gene_cat <- gene_cat[!duplicated(gene_cat)]
rp <- merge(rp, gene_cat, by.x = "geneID", by.y = "KNOWN_MUTATION")
rm(gene_cat)


# Rename categories 
rp[CATEGORY == "I", CATEGORY := "Pathogenic variants"]
rp[CATEGORY == "II", CATEGORY := "VUS"]


# Annotate with causal gene
rp$causal_gene <- NULL
rp[, causal_gene :=  sample_gene %in% sa$sample_gene]
rp[, cat3 := sample_gene %in% sa[CATEGORY == "III"]$sample_gene ]

rp[geneID == "STIM1", geneID := "STIM1*"]
rp[causal_gene == T & is.na(PROTEIN_ZSCORE), PROTEIN_ZSCORE := -Inf]
rp[causal_gene == T & is.na(RNA_ZSCORE), RNA_ZSCORE := -Inf]

rp_c <- rp[causal_gene == T & CATEGORY %in% c("Pathogenic variants" ,  "VUS")]
rp$geneID <- factor(rp$geneID, unique(rp_c[ order(PROTEIN_ZSCORE)]$geneID))


#####################################################################


#+ fig.width=10, fig.height=4
protrider <- ggplot(rp, aes(geneID,  PROTEIN_ZSCORE ))+
  geom_hline(yintercept = 0, color = "grey50") +
  geom_quasirandom(data = rp[causal_gene != T ] , aes(geneID,  PROTEIN_ZSCORE),color = "grey80" , size=0.5)+
  stat_summary(fun.data = mean_sdl, fun.args = list(mult=1), geom="errorbar", color="black", width=0.2) +
  stat_summary(fun = mean, geom="point", color="black") +
  ylab("PROTRIDER \n Protein zScore") + 
  geom_point( data = rp[causal_gene == T ], aes(geneID,  PROTEIN_ZSCORE), colour = "#FB9A99" , size = 2)+
  geom_point( data = rp[cat3 == T ], aes(geneID,  PROTEIN_ZSCORE ), size = 1, colour = "darkgreen", shape = 8)+
  theme_bw()+
  theme( legend.title = element_blank(),
         axis.title.x = element_blank(),
         plot.title = element_text(face="bold",  size=14, hjust = 0.5),
         axis.title.y = element_text(face="bold",  size=10, hjust = 0.5),
         axis.text.x = element_text(face="bold",  size=9, angle = 45, hjust = 1),
         axis.text.y = element_text(face="bold",  size=9),
         legend.position = "none")+
  facet_grid(~CATEGORY, scales = "free_x")


#+ fig.width=10, fig.height=4
outrider <- ggplot(rp, aes(geneID,  RNA_ZSCORE ))+
  geom_hline(yintercept = 0, color = "grey50") +
  geom_quasirandom(data = rp[causal_gene != T ] , aes(geneID,  RNA_ZSCORE),color = "grey80" , size=0.5)+
  stat_summary(fun.data = mean_sdl, fun.args = list(mult=1), geom="errorbar", color="black", width=0.2) +
  stat_summary(fun = mean, geom="point", color="black") +
  ylab("OUTRIDER \n RNA zScore") + 
  geom_point( data = rp[causal_gene == T ], aes(geneID,  RNA_ZSCORE), color = "#A6CEE3" , size = 2)+
  geom_point( data = rp[cat3 == T ], aes(geneID,  RNA_ZSCORE ), size = 1, colour = "darkgreen", shape = 8)+
  theme_bw()+
  theme( legend.title = element_blank(),
         axis.title.x = element_blank(),
         plot.title = element_text(face="bold",  size=14, hjust = 0.5),
         axis.title.y = element_text(face="bold",  size=10, hjust = 0.5),
         axis.text.x = element_text(face="bold",  size=9, angle = 45, hjust = 1),
         axis.text.y = element_text(face="bold",  size=9),
         legend.position = "none")+
  facet_grid(~CATEGORY, scales = "free_x")



#+ fig.width=10, fig.height=7
s_fig <- protrider / outrider
s_fig




pdf(snakemake@output$fig,
    width = 10, height =7,  useDingbats=FALSE )
print(s_fig) 
dev.off()


