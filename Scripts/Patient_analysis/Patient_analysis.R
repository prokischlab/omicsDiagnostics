#'---
#' title: Patient analysis 
#' author: smirnovd
#' wb:
#'  input:
#'  - sample_annotation: '`sm config["ANNOTATION"]`'
#'  - patient_omics: '`sm config["PROC_DATA"] + "/integration/patient_omics_full.RDS"`'
#'  - patient_hpo: '`sm config["RAW_HPO"]`'
#'  - patient_organ_syst: '`sm config["PROC_DATA"] + "/HPO/Patients_affected_organs.tsv"`'
#'  - complex_results: '`sm config["PROC_DATA"] + "/Complexes/Complex_outliers_LIMMA.rds"`'
#'  output:
#' output: 
#'   html_document:
#'    code_folding: hide
#'    code_download: TRUE
#'---

# Load config
source("src/config.R")

# Load plotting functions
source("src/functions/plots.R")

# Load sample annotation
# sa <- fread('/s/project/mitoMultiOmics/multiOMICs_integration/raw_data/proteomics_annotation.tsv')
sa <- fread(snakemake@input$sample_annotation)
sa <- sa[USE_FOR_PROTEOMICS_PAPER == T]
# Subset diagnosed cases and candidates
sa <- sa[ !is.na(KNOWN_MUTATION)]


# Read integrated omics file 
# rp <- readRDS("/s/project/mitoMultiOmics/multiOMICs_integration/processed_data/integration/patient_omics_full.RDS") %>% as.data.table()
rp <- readRDS(snakemake@input$patient_omics) %>% as.data.table()
# Subset diagnosed cases and candidates
rp <- rp[SAMPLE_ID %in% sa$SAMPLE_ID]
rpx <- rp
# Subset only outliers or causal/candidate genes
rp <- rp[ outlier_class != 'non_outlier' | causal_gene == T]



# complexes <- readRDS("/s/project/mitoMultiOmics/multiOMICs_integration/processed_data/Complexes/Complex_outliers_LIMMA.rds") %>% as.data.table()
complexes <- readRDS(snakemake@input$complex_results) %>% as.data.table()
# Subset diagnosed cases and candidates
complexes <- complexes[SAMPLE_ID %in% sa$SAMPLE_ID]



# Load patient's HPO 
pat_hpo <- fread(snakemake@input$patient_hpo)
# pat_hpo <- fread('/s/project/mitoMultiOmics/multiOMICs_integration/raw_data/Patient_HPO_phenotypes.tsv')

# Subset diagnosed cases and candidates
pat_hpo <- pat_hpo[SAMPLE_ID %in% sa$SAMPLE_ID, c("SAMPLE_ID", "HPO_term")]
pat_hpo <-pat_hpo[!duplicated(pat_hpo ),]

# Load patient's affected organ systems
patients <- fread(snakemake@input$patient_organ_syst)
# patients <- fread('/s/project/mitoMultiOmics/multiOMICs_integration/processed_data/HPO/Patients_affected_organs.tsv')
# Subset diagnosed cases and candidates
patients <- patients[SAMPLE_ID %in% sa$SAMPLE_ID]

# i <- "OM94976" # for tests

#' ## Category I - Known pathogenic variants
#+ fig.width=18, fig.height=12
for (i in unique(sa[ CATEGORY == "I" ]$SAMPLE_ID)){
  print(paste("Sample ID:", i, "Causal gene:", sa[SAMPLE_ID == i]$KNOWN_MUTATION ) )
  
  plot_org <- gganatogram(data=  patients[SAMPLE_ID == i], 
                          fillOutline='white', 
                          organism='human', 
                          sex= unique(sa[SAMPLE_ID == i]$gender), 
                          fill="colour")+ 
    theme_void() + ggtitle(paste0("individual #", i )) +
    theme(plot.title = element_text(margin = NULL,face="bold", size = 9, hjust = 0.5),
          plot.margin = margin(0, 0, 0, 0, "cm"))
  
  phenotypes <- pat_hpo[SAMPLE_ID == i,  c( "HPO_term") ]
  phenotypes <- gridExtra::tableGrob(phenotypes)
  plot_org <-  plot_org  | phenotypes 
  
  if (nrow(rp[SAMPLE_ID == i ]) == 0){
    print("No RNA-seq avaliable")
    # fig <- protein_volcanoZ(rpx, sample = i)
    plot_sample <- plot_org

  }else{
    plot_sample <- plot_patient_all(rp, sample = i) 
    plot_sample <- plot_org | plot_sample
  }
  
  volcano_rna <- rna_volcanoF(rpx, i) + theme( legend.position = "none") + ggtitle("Aberrant RNA expression")
  volcano_protein <- protein_volcanoF(rpx, i) + theme( legend.position = "none")  + ggtitle("Aberrant protein expression")
  volcano_complex <- complex_volcanoF(complexes, i)  + ggtitle("Aberrant protein complex expression")
  
  volcanos <- volcano_rna | volcano_protein | volcano_complex
  
  plot_sample <- plot_sample / volcanos
  print(plot_sample)
}




#' ## Category II - VUS

#' ### Category IIa - VUS validated
#+ fig.width=18, fig.height=12
for (i in unique(sa[ CATEGORY == "IIa" ]$SAMPLE_ID)){
  print(paste("Sample ID:", i, "Causal gene:", sa[SAMPLE_ID == i]$KNOWN_MUTATION ) )
  
  plot_org <- gganatogram(data=  patients[SAMPLE_ID == i], 
                          fillOutline='white', 
                          organism='human', 
                          sex= unique(sa[SAMPLE_ID == i]$gender), 
                          fill="colour")+ 
    theme_void() + ggtitle(paste0("individual #", i )) +
    theme(plot.title = element_text(margin = NULL,face="bold", size = 9, hjust = 0.5),
          plot.margin = margin(0, 0, 0, 0, "cm"))
  
  phenotypes <- pat_hpo[SAMPLE_ID == i,  c( "HPO_term") ]
  phenotypes <- gridExtra::tableGrob(phenotypes)
  plot_org <-  plot_org  | phenotypes 
  
  if (nrow(rp[SAMPLE_ID == i ]) == 0){
    print("No outliers detected")
    plot_sample <- plot_org
  }else{
    plot_sample <- plot_patient_all(rp, sample = i) 
    plot_sample <- plot_org | plot_sample
  }
  
  
  volcano_rna <- rna_volcanoF(rpx, i) + theme( legend.position = "none") + ggtitle("Aberrant RNA expression")
  volcano_protein <- protein_volcanoF(rpx, i) + theme( legend.position = "none")  + ggtitle("Aberrant protein expression")
  volcano_complex <- complex_volcanoF(complexes, i)  + ggtitle("Aberrant protein complex expression")
  
  volcanos <- volcano_rna | volcano_protein | volcano_complex
  
  plot_sample <- plot_sample / volcanos
  
  print(plot_sample)
}


#' ### Category IIb - VUS not validated
#+ fig.width=18, fig.height=12
for (i in unique(sa[ CATEGORY == "IIb" ]$SAMPLE_ID)){
  print(paste("Sample ID:", i, "Causal gene:", sa[SAMPLE_ID == i]$KNOWN_MUTATION ) )
  
  plot_org <- gganatogram(data=  patients[SAMPLE_ID == i], 
                          fillOutline='white', 
                          organism='human', 
                          sex= unique(sa[SAMPLE_ID == i]$gender), 
                          fill="colour")+ 
    theme_void() + ggtitle(paste0("individual #", i )) +
    theme(plot.title = element_text(margin = NULL,face="bold", size = 9, hjust = 0.5),
          plot.margin = margin(0, 0, 0, 0, "cm"))
  
  phenotypes <- pat_hpo[SAMPLE_ID == i,  c( "HPO_term") ]
  phenotypes <- gridExtra::tableGrob(phenotypes)
  plot_org <-  plot_org  | phenotypes 
  
  if (nrow(rp[SAMPLE_ID == i ]) == 0){
    print("No outliers detected")
    plot_sample <- plot_org
    
  }else{
    plot_sample <- plot_patient_all(rp, sample = i) 
    plot_sample <- plot_org | plot_sample
  }
  
  
  volcano_rna <- rna_volcanoF(rpx, i) + theme( legend.position = "none") + ggtitle("Aberrant RNA expression")
  volcano_protein <- protein_volcanoF(rpx, i) + theme( legend.position = "none")  + ggtitle("Aberrant protein expression")
  volcano_complex <- complex_volcanoF(complexes, i)  + ggtitle("Aberrant protein complex expression")
  
  volcanos <- volcano_rna | volcano_protein | volcano_complex
  
  plot_sample <- plot_sample / volcanos
  
  print(plot_sample)
}




#' ### Category IIc - VUS not detected
#+ fig.width=18, fig.height=12
for (i in unique(sa[ CATEGORY == "IIc" ]$SAMPLE_ID)){
  print(paste("Sample ID:", i, "Causal gene:", sa[SAMPLE_ID == i]$KNOWN_MUTATION ) )
  
  plot_org <- gganatogram(data=  patients[SAMPLE_ID == i], 
                          fillOutline='white', 
                          organism='human', 
                          sex= unique(sa[SAMPLE_ID == i]$gender), 
                          fill="colour")+ 
    theme_void() + ggtitle(paste0("individual #", i )) +
    theme(plot.title = element_text(margin = NULL,face="bold", size = 9, hjust = 0.5),
          plot.margin = margin(0, 0, 0, 0, "cm"))
  
  phenotypes <- pat_hpo[SAMPLE_ID == i,  c( "HPO_term") ]
  phenotypes <- gridExtra::tableGrob(phenotypes)
  plot_org <-  plot_org  | phenotypes 
  
  if (nrow(rp[SAMPLE_ID == i ]) == 0){
    print("No outliers detected")
    plot_sample <- plot_org
    
  }else{
    plot_sample <- plot_patient_all(rp, sample = i) 
    plot_sample <- plot_org | plot_sample
  }
  
  
  volcano_rna <- rna_volcanoF(rpx, i) + theme( legend.position = "none") + ggtitle("Aberrant RNA expression")
  volcano_protein <- protein_volcanoF(rpx, i) + theme( legend.position = "none")  + ggtitle("Aberrant protein expression")
  volcano_complex <- complex_volcanoF(complexes, i)  + ggtitle("Aberrant protein complex expression")
  
  volcanos <- volcano_rna | volcano_protein | volcano_complex
  
  plot_sample <- plot_sample / volcanos
  
  print(plot_sample)
}


#' ## Category III New discoveries
#+ fig.width=18, fig.height=12
for (i in unique(sa[ CATEGORY == "III" ]$SAMPLE_ID)){
  print(paste("Sample ID:", i, "Causal gene:", sa[SAMPLE_ID == i]$KNOWN_MUTATION ) )
  
  plot_org <- gganatogram(data=  patients[SAMPLE_ID == i], 
                          fillOutline='white', 
                          organism='human', 
                          sex= unique(sa[SAMPLE_ID == i]$gender), 
                          fill="colour")+ 
    theme_void() + ggtitle(paste0("individual #", i )) +
    theme(plot.title = element_text(margin = NULL,face="bold", size = 9, hjust = 0.5),
          plot.margin = margin(0, 0, 0, 0, "cm"))
  
  phenotypes <- pat_hpo[SAMPLE_ID == i,  c( "HPO_term") ]
  phenotypes <- gridExtra::tableGrob(phenotypes)
  plot_org <-  plot_org  | phenotypes 
  
  if (nrow(rp[SAMPLE_ID == i ]) == 0){
    print("No outliers detected")
    plot_sample <- plot_org
    
  }else{
    plot_sample <- plot_patient_all(rp, sample = i) 
    plot_sample <- plot_org | plot_sample
  }
  
  
  volcano_rna <- rna_volcanoF(rpx, i) + theme( legend.position = "none") + ggtitle("Aberrant RNA expression")
  volcano_protein <- protein_volcanoF(rpx, i) + theme( legend.position = "none")  + ggtitle("Aberrant protein expression")
  volcano_complex <- complex_volcanoF(complexes, i)  + ggtitle("Aberrant protein complex expression")
  
  volcanos <- volcano_rna | volcano_protein | volcano_complex
  
  plot_sample <- plot_sample / volcanos
  
  print(plot_sample)
}




