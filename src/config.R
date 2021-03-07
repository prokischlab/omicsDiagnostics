suppressPackageStartupMessages({
  stringsAsFactors=F
  library(vctrs)
  library(data.table)
  library(readxl)
  library(readr)
 
  library(plyr)
  library(dplyr)
  library(devtools)
  library(tidyr)
  library(tidyverse)
  library(stringr)
  
  library(ggplot2)
  library(ggrepel)
  library(plotly)
  library(ggforce)
  library(gplots)
  library(ggthemes)
  library(pheatmap)
  library(ggbeeswarm)
  library(patchwork)
  library(cowplot)
  library(gganatogram) 
  library(ggpubr)

  library(reshape2)
  library(limma)
  library(MASS)
  library(biomaRt)
  library(magrittr)
  library(seqinr)

  library(LSD)
  library(knitr)
  library(markdown)

  library(SummarizedExperiment)
  library(OUTRIDER)
  library(clusterProfiler)

 
  library(ontologyIndex)
  library(PCAN)
  library(ontologySimilarity)
  
})


# Outlier threshold definition 
PADJ_THRESHOLD =  0.1
# Outlier threshold definition for rna-and-protein outliers
ZSCORE_THRESHOLD = 3
LOG2FC_THRESHOLD = 1

MIN_SAMPLE_NUM_PROT <- 0.52 #0.76 # % of samples protein not detected in


# Define colours for plotting 

outlier_colors <- c(  "non_outlier" = "gray80",  
                      "protein-only" = "#FB9A99",  
                      "RNA-only" = "#A6CEE3",       
                      "RNA-and-protein" = "#B2DF8A",  
                      "no rare variant" = "grey20", 
                      "1 rare variant" = "black",
                      "RNA_underexpression" = "#A6CEE3", 
                      "Protein_underexpression" = "#FB9A99", 
                      "RNA_Protein_underexpression" = "#B2DF8A",
                      "RNA_overexpression" = "blue", 
                      "Protein_overexpression" = "red", 
                      "RNA_Protein_overexpression" = "green")


variant_colors <- c("non_coding" = "gray81",  
                    "splice" = "darkolivegreen2" ,  
                    "stop" =  "aquamarine3",  
                    "coding" = "salmon1",   
                    "frameshift" = "lightskyblue3",  
                    "synonymous" = "moccasin",
                    "no rare variant" = "white",
                    "rare" = "white", 
                    "rare " = "white", 
                    "potential_biallelic" = "gray70",
                    "potential biallelic SemSim > 2" = "gray50",
                    "potential biallelic + Semantic sim > 2" = "gray50")

variant_colors2 <- c("no rare variant" = "white",
                    "SemSim > 2"= "gray90",
                    "rare" = "gray70", 
                    "rare SemSim > 2" = "gray50",
                    "potential_biallelic" = "gray30",
                    "potential biallelic SemSim > 2" = "gray10") 


gene_shapes <- c("no data" = 20, 
                 "no rare variant" = 20,  
                 "1 rare variant" = 17, 
                 "rare pot. biallelic variants" = 15) 

gene_size <- c("no data" = 2.5, 
               "no rare variant" = 2.5, 
                "1 rare variant" = 3.1,
                "rare pot. biallelic variants" = 4.5) 

text_color <- c("no data" = "grey20", 
                "no rare variant" = "grey20", 
                "1 rare variant" = "black", 
                "rare pot. biallelic variants" = "black")










