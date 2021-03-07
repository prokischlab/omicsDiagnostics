#'---
#' title: Supplementary Fig 3 
#' author: smirnovd
#' wb:
#'  input: 
#'  - config: 'src/config.R'
#'  - patient_omics: '`sm config["PROC_DATA"] + "/integration/patient_omics_full.RDS"`'
#'  - sample_annotation: '`sm config["ANNOTATION"]`'
#'  output:
#'  - fig: '`sm config["FIGURE_DIR"] + "/Supplementary_figures/S_Fig3.pdf"`'
#' output: 
#'   html_document:
#'    code_folding: hide
#'    code_download: TRUE
#'---


# Load config
source(snakemake@input$config)

# Load plotting functions
source("src/functions/plots.R")


# Load sample annotation
# sa <- fread('/s/project/mitoMultiOmics/multiOMICs_integration/raw_data/proteomics_annotation.tsv')
sa <- fread(snakemake@input$sample_annotation)
sa <- sa[USE_FOR_PROTEOMICS_PAPER == T]

# Subset diagnosed cases from discovery cohort
sa <- sa[ !is.na(KNOWN_MUTATION) & CATEGORY == "III"]

# Remove cases appearing in the main figures
sa <- sa[ !(KNOWN_MUTATION %in% c("LIG3", "MRPL38", "MORC2", "DARS2", "MRPS25", "EPG5") )]

# Read integrated omics file 
# rp <- readRDS("/s/project/mitoMultiOmics/multiOMICs_integration/processed_data/integration/patient_omics.RDS") %>% as.data.table()
rp <- readRDS(snakemake@input$patient_omics) %>% as.data.table()
# Subset diagnosed cases and candidates
rp <- rp[SAMPLE_ID %in% sa$SAMPLE_ID]

# Subset outliers only 
rp <- rp[ outlier_class != 'non_outlier' | causal_gene == T]




samples <- unique(sa$SAMPLE_ID)

# #+ fig.width=5, fig.height=5
# for (i in samples){
#   plot_sample <- plot_patient_main(rp, sample = i) + 
#     ggtitle(paste("Sample ID:", i, "\nCasual gene:", sa[SAMPLE_ID == i]$KNOWN_MUTATION ) )
#   print(plot_sample)
# }



 
plot_OM87369 <- plot_patient_main(rp, sample = "OM87369") +  ggtitle(paste("Sample ID:", "OM87369", "\nMolecular diagnosis:", sa[SAMPLE_ID == "OM87369"]$KNOWN_MUTATION ) ) #+ theme(legend.title = element_blank(), legend.position = "none") 
plot_OM34814 <- plot_patient_main(rp, sample = "OM34814") +  ggtitle(paste("Sample ID:", "OM34814", "\nMolecular diagnosis:", sa[SAMPLE_ID == "OM34814"]$KNOWN_MUTATION ) ) #+ theme(legend.title = element_blank(), legend.position = "none") 
plot_OM38813 <- plot_patient_main(rp, sample = "OM38813") +  ggtitle(paste("Sample ID:", "OM38813", "\nMolecular diagnosis:", sa[SAMPLE_ID == "OM38813"]$KNOWN_MUTATION ) ) #+ theme(legend.title = element_blank(), legend.position = "none") 
plot_OM21111 <- plot_patient_main(rp, sample = "OM21111") +  ggtitle(paste("Sample ID:", "OM21111", "\nMolecular diagnosis:", sa[SAMPLE_ID == "OM21111"]$KNOWN_MUTATION ) ) #+ theme(legend.title = element_blank(), legend.position = "none") 
plot_OM56706 <- plot_patient_main(rp, sample = "OM56706") +  ggtitle(paste("Sample ID:", "OM56706", "\nMolecular diagnosis:", sa[SAMPLE_ID == "OM56706"]$KNOWN_MUTATION ) ) #+ theme(legend.title = element_blank(), legend.position = "none") 
plot_OM03592 <- plot_patient_main(rp, sample = "OM03592") +  ggtitle(paste("Sample ID:", "OM03592", "\nMolecular diagnosis:", sa[SAMPLE_ID == "OM03592"]$KNOWN_MUTATION ) )

s_fig <- plot_OM87369 + plot_OM34814 + plot_OM38813 +  plot_OM21111 + plot_OM56706  + plot_OM03592  + 
  plot_layout(ncol = 2) # + plot_annotation(tag_levels = 'a') & theme(plot.tag = element_text(size = 14, face = "bold") )

# #+ fig.width=13, fig.height=16
s_fig


pdf( snakemake@output$fig, # "/s/project/mitoMultiOmics/multiOMICs_integration/Figures/Supplementary_figures/S_Fig3.pdf",   # 
    useDingbats=FALSE, width = 13, height =16,  ) # paper = "a4"
s_fig
dev.off()








