#########################
### plotting functions###
#########################

suppressPackageStartupMessages({
  stringsAsFactors=F
  require(ggplot2)
  require(patchwork)
  require(ggrepel)
  require(gganatogram)
  require(gridExtra)
})


plot_patient_all <- function(data, sample = NULL) {
  require(patchwork)
  require(data.table)
  
  samp <- data[ SAMPLE_ID == sample] %>% as.data.table() #subset for patient
  range <- c(samp$RNA_ZSCORE, samp$PROTEIN_ZSCORE)
  
  samp[is.na(RNA_ZSCORE)]$RNA_ZSCORE <- 0
  samp[is.na(PROTEIN_ZSCORE)]$PROTEIN_ZSCORE <- 0
  
  samp <- samp[ !duplicated(samp[, c("SAMPLE_ID", "geneID")]) ]
  samp[ is.na(Semantic_sim), Semantic_sim:= 0]
  
  main <- ggplot(samp[gene_detected == "RNA and protein detected" ], aes(RNA_ZSCORE, PROTEIN_ZSCORE)) +
    geom_vline(xintercept = 0, color = "grey50") +
    geom_hline(yintercept = 0, color = "grey50") +
    
    geom_point(aes(color= outlier_class, shape= gene_class, size = Semantic_sim )) + 
    
    scale_color_manual(breaks = c("non_outlier","RNA-only", "protein-only", "RNA-and-protein"),
                       values = outlier_colors)+
    
    scale_shape_manual(values = gene_shapes) +
    
    
    geom_text_repel(data= samp[gene_detected == "RNA and protein detected" & (gene_class == "no rare variant" | gene_class == "no data") ], 
                    mapping=aes(RNA_ZSCORE, PROTEIN_ZSCORE, label= geneID), colour= "grey20", 
                    size= 2.8, show.legend = F, min.segment.length = Inf)+  # , min.segment.length = Inf
    
    geom_text_repel(data= samp[gene_detected == "RNA and protein detected" & gene_class == "1 rare variant" ], 
                    mapping=aes(RNA_ZSCORE, PROTEIN_ZSCORE, label= geneID), colour= "black", 
                    size= 3.3, show.legend = F)+ 
    
    geom_label_repel(data= samp[gene_detected == "RNA and protein detected" & gene_class == "rare pot. biallelic variants"], 
                     mapping=aes(RNA_ZSCORE, PROTEIN_ZSCORE, label= geneID ), 
                     box.padding = unit(0.3, "lines"), colour= "black", 
                     size= 4, show.legend = F )+ 
    
    coord_fixed( xlim = c( -max( abs(range) , na.rm = T), max( abs(range) , na.rm = T) ),  
                 ylim = c( -max( abs(range) , na.rm = T),  max( abs(range) , na.rm = T) )) +  
    theme_bw()+
    labs(color='Outlier class',
         size='Semantic similarity',
         shape='Variant information') +
    scale_size(range = c(2.5,6))+
    theme( plot.title = element_text(hjust = 0.5),
           axis.title.x= element_blank(),
           axis.title.y= element_blank(),
           #legend.title = element_blank(),
           axis.text.x = element_text(face="bold",  size=10),
           axis.text.y = element_text(face="bold",  size=10),
           plot.margin = margin(0, 0, 0, 0, "cm") ) 
  
  
  rna <- ggplot(samp[gene_detected == "no Protein" ], aes(RNA_ZSCORE, PROTEIN_ZSCORE)) +
    geom_vline(xintercept = 0, color = "grey50") +
    geom_point(aes(color= outlier_class, shape= gene_class, size = Semantic_sim)) + 
    
    xlab("RNA zScore") + 
    ylab("Not detected\n by proteomics") +
    
    scale_color_manual(breaks = c("non_outlier","RNA-only","protein-only", "RNA-and-protein"),
                       values = outlier_colors)+
    scale_shape_manual(values = gene_shapes) +
    
    
    geom_text_repel(data= samp[gene_detected == "no Protein" & gene_class == "no rare variant" ], 
                    mapping=aes(RNA_ZSCORE, PROTEIN_ZSCORE, label= geneID), colour= "grey20", size= 2.8, show.legend = F)+ 
    
    geom_text_repel(data= samp[gene_detected == "no Protein" & gene_class == "1 rare variant" ], 
                    mapping=aes(RNA_ZSCORE, PROTEIN_ZSCORE, label= geneID), colour= "black", size= 3.3, show.legend = F)+ 
    
    geom_label_repel(data= samp[gene_detected == "no Protein" & gene_class == "rare pot. biallelic variants"], 
                     mapping=aes(RNA_ZSCORE, PROTEIN_ZSCORE, label= geneID ), 
                     box.padding = unit(0.3, "lines"), colour= "black", size= 4, show.legend = F )+ 
    
    
    coord_fixed( xlim = c( -max( abs(range) , na.rm = T), max( abs(range) , na.rm = T)  ), 
                 ylim = c( -0.7,  0.7  )) +  
    theme_linedraw() +
    scale_size(range = c(2.5,6))+
    theme(axis.title.y= element_text( size=6, margin = NULL,face="bold"),
          axis.title.x= element_text( size=10, margin = NULL,face="bold"),
          axis.text.x = element_text(face="bold",  size=10),
          axis.text.y=element_blank(), 
          legend.title = element_blank(),
          axis.ticks.y=element_blank(),
          legend.position = "none",
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          plot.margin = margin(0, 0, 0, 0, "cm"))
  
  
  protein <- ggplot(samp[gene_detected == "no RNA" ], aes(RNA_ZSCORE, PROTEIN_ZSCORE)) +
    geom_hline(yintercept = 0, color = "grey50") +
    geom_point(aes(color= outlier_class, shape= gene_class, size = Semantic_sim)) + 
    
    xlab("Not detected\n by RNA-seq") + 
    ylab("Protein zScore") +
    
    scale_color_manual(breaks = c("non_outlier","RNA-only","protein-only", "RNA-and-protein"),
                       values = outlier_colors)+
    scale_shape_manual(values = gene_shapes) +
    
    geom_text_repel(data= samp[gene_detected == "no RNA" & gene_class == "no rare variant" ], 
                    mapping=aes(RNA_ZSCORE, PROTEIN_ZSCORE, label= geneID), colour= "grey20", size= 2.8, show.legend = F)+ 
    
    geom_text_repel(data= samp[gene_detected == "no RNA" & gene_class == "1 rare variant" ], 
                    mapping=aes(RNA_ZSCORE, PROTEIN_ZSCORE, label= geneID), colour= "black", size= 3.3, show.legend = F)+ 
    
    geom_label_repel(data= samp[gene_detected == "no RNA" & gene_class == "rare pot. biallelic variants"], 
                     mapping=aes(RNA_ZSCORE, PROTEIN_ZSCORE, label= geneID ), 
                     box.padding = unit(0.3, "lines"), colour= "black", size= 4, show.legend = F )+ 
    
    
    coord_fixed( ylim = c( -max( abs(range) , na.rm = T), max( abs(range) , na.rm = T) ), 
                 xlim = c( -0.7,  0.7  )) + 
    
    theme_linedraw() +
    scale_size(range = c(2.5,6))+
    theme(axis.title.x= element_text( size=6, margin = NULL,face="bold" ), 
          axis.title.y= element_text( size=10, margin = NULL,face="bold"),
          axis.text.y = element_text(face="bold",  size=10),
          axis.text.x = element_blank(), 
          legend.title = element_blank(),
          axis.ticks.x = element_blank(),
          legend.position = "none",
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          plot.margin = margin(0, 0, 0, 0, "cm"))
  
  layout <- '
AB
#C
'
  figure <- wrap_plots(A = protein, B = main, C = rna, design = layout)
  return(figure)
}


plot_patient_main <- function(data, sample = NULL) {
  require(data.table)
  samp <- data[ SAMPLE_ID == sample & gene_detected == "RNA and protein detected"] %>% as.data.table() #subset for patient
  samp <- samp[!is.na(RNA_ZSCORE)  & !is.na(PROTEIN_ZSCORE)] #remove missing rna and proteins
  
  range <- c(samp$RNA_ZSCORE, samp$PROTEIN_ZSCORE)
  samp <- samp[ !duplicated(samp[, c("SAMPLE_ID", "geneID")]) ]
  samp[ is.na(Semantic_sim), Semantic_sim:= 0]
  
  
  figure <- ggplot(samp, aes(RNA_ZSCORE, PROTEIN_ZSCORE)) +
    geom_vline(xintercept = 0, color = "grey50") +
    geom_hline(yintercept = 0, color = "grey50") +
    
    geom_point(aes(color= outlier_class, shape= gene_class, size = Semantic_sim )) + 
    
    xlab("RNA zScore") + 
    ylab("Protein zScore") +
    scale_color_manual(breaks = c("non_outlier","RNA-only","protein-only", "RNA-and-protein"),
                       values = outlier_colors)+
    
    scale_shape_manual(values = gene_shapes) +
    
    
    
    geom_text_repel(data= samp[gene_class == "no rare variant" | gene_class == "no data" ], 
                    mapping=aes(RNA_ZSCORE, PROTEIN_ZSCORE, label= geneID), colour= "grey20", 
                    size= 2.8, show.legend = F, min.segment.length = Inf)+  # , min.segment.length = Inf 
    
    geom_text_repel(data= samp[gene_class == "1 rare variant" ], 
                    mapping=aes(RNA_ZSCORE, PROTEIN_ZSCORE, label= geneID), colour= "black", 
                    size= 3.3, show.legend = F)+ 
    
    geom_label_repel(data= samp[gene_class == "rare pot. biallelic variants"], 
                     mapping=aes(RNA_ZSCORE, PROTEIN_ZSCORE, label= geneID ), 
                     box.padding = unit(0.3, "lines"), colour= "black", 
                     size= 4, show.legend = F )+ 
    
    
    coord_fixed( xlim = c( -max( abs(range) , na.rm = T), max( abs(range) , na.rm = T)  ),  
                 ylim = c( -max( abs(range) , na.rm = T),  max( abs(range) , na.rm = T)  )) +  
    theme_bw()+
    labs(color='Outlier class',
         size='Semantic similarity',
         shape='Variant information') +
    scale_size(range = c(2.5,6))+
    theme(plot.title = element_text(hjust = 0.5),
          axis.title.x= element_text( size=10, margin = NULL,face="bold"),
          axis.title.y= element_text( size=10, margin = NULL,face="bold"),
          axis.text.x = element_text(face="bold",  size=10),
          axis.text.y = element_text(face="bold",  size=10),
          plot.margin = margin(0, 0, 0, 0, "cm") ) 
  
  return(figure)
}



###########################################################
###########################################################

protein_volcanoZ <- function(data, sample = NULL) {
  require(data.table)
  samp <- data[ SAMPLE_ID == sample] %>% as.data.table() #subset for patient
  samp <- samp[!is.na(PROTEIN_ZSCORE)]  #remove missing
  samp <- samp[ !duplicated(samp[, c("SAMPLE_ID", "geneID")]) ]
  samp[, outlier_class := "non_outlier"]
  samp[PROTEIN_outlier == T, outlier_class := "protein-only"]
  samp[ is.na(Semantic_sim), Semantic_sim:= 0]
  
  figure <- ggplot(samp, aes(PROTEIN_ZSCORE, -log10(PROTEIN_PVALUE) )) +
    geom_vline(xintercept = 0, color = "grey50") +
    
    geom_point(aes(color= outlier_class, shape= gene_class, size = Semantic_sim )) + 
    
    xlab("Protein zScore") + 
    ylab("-log10(P-value)") +
    scale_color_manual(breaks = c("non_outlier","RNA-only","protein-only", "RNA-and-protein"),
                       values = outlier_colors)+
    
    scale_shape_manual(values = gene_shapes) +
    
    
    
    geom_text_repel(data= samp[outlier_class != "non_outlier" & (gene_class == "no rare variant" | gene_class == "no data") ], 
                    mapping=aes(PROTEIN_ZSCORE, -log10(PROTEIN_PVALUE), label= geneID), colour= "grey20", 
                    size= 2.8, show.legend = F, min.segment.length = Inf)+  # , min.segment.length = Inf
    
    geom_text_repel(data= samp[outlier_class != "non_outlier" & gene_class == "1 rare variant" ], 
                    mapping=aes(PROTEIN_ZSCORE, -log10(PROTEIN_PVALUE), label= geneID), colour= "black", 
                    size= 3.3, show.legend = F)+ 
    
    geom_label_repel(data= samp[outlier_class != "non_outlier" & gene_class == "rare pot. biallelic variants"], 
                     mapping=aes(PROTEIN_ZSCORE, -log10(PROTEIN_PVALUE), label= geneID ), 
                     box.padding = unit(0.3, "lines"), colour= "black", 
                     size= 4, show.legend = F )+ 
    
    
    theme_bw()+
    scale_size(range = c(2.5,6))+
    theme(plot.title = element_text(hjust = 0.5),
          axis.title.x= element_text( size=10, margin = NULL,face="bold"),
          axis.title.y= element_text( size=10, margin = NULL,face="bold"),
          axis.text.x = element_text(face="bold",  size=10),
          axis.text.y = element_text(face="bold",  size=10),
          legend.title = element_blank(),
          plot.margin = margin(0, 0, 0, 0, "cm") ) 
  
  return(figure)
}


rna_volcanoZ <- function(data, rna = NULL) {
  require(data.table)
  samp <- data[ SAMPLE_ID == rna] %>% as.data.table() #subset for patient
  samp <- samp[!is.na(RNA_ZSCORE)]  #remove missing 
  samp <- samp[ !duplicated(samp[, c("SAMPLE_ID", "geneID")]) ]
  samp[, outlier_class := "non_outlier"]
  samp[RNA_outlier == T, outlier_class := "RNA-only"]
  samp[ is.na(Semantic_sim), Semantic_sim:= 0]
  
  figure <- ggplot(samp, aes(RNA_ZSCORE, -log10(RNA_PVALUE) )) +
    geom_vline(xintercept = 0, color = "grey50") +
    geom_point(aes(color= outlier_class, shape= gene_class, size = Semantic_sim )) + 
    
    xlab("RNA zScore") + 
    ylab("-log10(P-value)") +
    scale_color_manual(breaks = c("non_outlier","RNA-only","protein-only", "RNA-and-protein"),
                       values = outlier_colors)+
    
    scale_shape_manual(values = gene_shapes) +
    
    
    
    geom_text_repel(data= samp[outlier_class != "non_outlier" & (gene_class == "no rare variant" | gene_class == "no data")  ], 
                    mapping=aes(RNA_ZSCORE, -log10(RNA_PVALUE), label= geneID), colour= "grey20", 
                    size= 2.8, show.legend = F, min.segment.length = Inf)+  # , min.segment.length = Inf
    
    geom_text_repel(data= samp[outlier_class != "non_outlier" & gene_class == "1 rare variant" ], 
                    mapping=aes(RNA_ZSCORE, -log10(RNA_PVALUE), label= geneID), colour= "black", 
                    size= 3.3, show.legend = F)+ 
    
    geom_label_repel(data= samp[outlier_class != "non_outlier" & gene_class == "rare pot. biallelic variants"], 
                     mapping=aes(RNA_ZSCORE, -log10(RNA_PVALUE), label= geneID ), 
                     box.padding = unit(0.3, "lines"), colour= "black", 
                     size= 4, show.legend = F )+ 
    
    
    theme_bw()+
    scale_size(range = c(2.5,6))+
    theme(plot.title = element_text(hjust = 0.5),
          axis.title.x= element_text( size=10, margin = NULL,face="bold"),
          axis.title.y= element_text( size=10, margin = NULL,face="bold"),
          axis.text.x = element_text(face="bold",  size=10),
          axis.text.y = element_text(face="bold",  size=10),
          legend.title = element_blank(),
          plot.margin = margin(0, 0, 0, 0, "cm") ) 
  
  return(figure)
}



###########################################################
###########################################################

protein_volcanoF <- function(data, sample = NULL) {
  require(data.table)
  samp <- data[ SAMPLE_ID == sample] %>% as.data.table() #subset for patient
  samp <- samp[!is.na(PROTEIN_LOG2FC)]  #remove missing
  samp <- samp[ !duplicated(samp[, c("SAMPLE_ID", "geneID")]) ]
  samp[, outlier_class := "non_outlier"]
  samp[PROTEIN_outlier == T, outlier_class := "protein-only"]
  samp[ is.na(Semantic_sim), Semantic_sim:= 0]
  
  figure <- ggplot(samp, aes(PROTEIN_LOG2FC, -log10(PROTEIN_PVALUE) )) +
    geom_vline(xintercept = 0, color = "grey50") +
    
    geom_point(aes(color= outlier_class, shape= gene_class, size = Semantic_sim )) + 
    
    xlab("Protein log2(Fold change)") + 
    ylab("-log10(P-value)") +
    scale_color_manual(breaks = c("non_outlier","RNA-only","protein-only", "RNA-and-protein"),
                       values = outlier_colors)+
    
    scale_shape_manual(values = gene_shapes) +
    
    
    
    geom_text_repel(data= samp[outlier_class != "non_outlier" & (gene_class == "no rare variant" | gene_class == "no data")  ], 
                    mapping=aes(PROTEIN_LOG2FC, -log10(PROTEIN_PVALUE), label= geneID), colour= "grey20", 
                    size= 2.8, show.legend = F, min.segment.length = Inf)+  # , min.segment.length = Inf
    
    geom_text_repel(data= samp[outlier_class != "non_outlier" & gene_class == "1 rare variant" ], 
                    mapping=aes(PROTEIN_LOG2FC, -log10(PROTEIN_PVALUE), label= geneID), colour= "black", 
                    size= 3.3, show.legend = F)+ 
    
    geom_label_repel(data= samp[ (outlier_class != "non_outlier" & gene_class == "rare pot. biallelic variants") | causal_gene == T ] , 
                     mapping=aes(PROTEIN_LOG2FC, -log10(PROTEIN_PVALUE), label= geneID ), 
                     box.padding = unit(0.3, "lines"), colour= "black", 
                     size= 4, show.legend = F )+ 
    
    
    theme_bw()+
    scale_size(range = c(2.5,6))+
    theme(plot.title = element_text(hjust = 0.5),
          axis.title.x= element_text( size=10, margin = NULL,face="bold"),
          axis.title.y= element_text( size=10, margin = NULL,face="bold"),
          axis.text.x = element_text(face="bold",  size=10),
          axis.text.y = element_text(face="bold",  size=10),
          legend.title = element_blank()) 
  
  return(figure)
}


rna_volcanoF <- function(data, rna = NULL) {
  require(data.table)
  samp <- data[ SAMPLE_ID == rna] %>% as.data.table() #subset for patient
  samp <- samp[!is.na(RNA_LOG2FC)]  #remove missing 
  samp <- samp[ !duplicated(samp[, c("SAMPLE_ID", "geneID")]) ]
  samp[, outlier_class := "non_outlier"]
  samp[RNA_outlier == T, outlier_class := "RNA-only"]
  samp[ is.na(Semantic_sim), Semantic_sim:= 0]
  
  figure <- ggplot(samp, aes(RNA_LOG2FC, -log10(RNA_PVALUE) )) +
    geom_vline(xintercept = 0, color = "grey50") +
    geom_point(aes(color= outlier_class, shape= gene_class, size = Semantic_sim )) + 
    
    xlab("RNA log2(Fold change)") + 
    ylab("-log10(P-value)") +
    scale_color_manual(breaks = c("non_outlier","RNA-only","protein-only", "RNA-and-protein"),
                       values = outlier_colors)+
    
    scale_shape_manual(values = gene_shapes) +
    
    
    
    geom_text_repel(data= samp[outlier_class != "non_outlier" & (gene_class == "no rare variant" | gene_class == "no data") ], 
                    mapping=aes(RNA_LOG2FC, -log10(RNA_PVALUE), label= geneID), colour= "grey20", 
                    size= 2.8, show.legend = F, min.segment.length = Inf)+  # , min.segment.length = Inf
    
    geom_text_repel(data= samp[outlier_class != "non_outlier" & gene_class == "1 rare variant" ], 
                    mapping=aes(RNA_LOG2FC, -log10(RNA_PVALUE), label= geneID), colour= "black", 
                    size= 3.4, show.legend = F)+ 
    
    geom_label_repel(data= samp[ (outlier_class != "non_outlier" & gene_class == "rare pot. biallelic variants") | causal_gene == T], 
                     mapping=aes(RNA_LOG2FC, -log10(RNA_PVALUE), label= geneID ), 
                     box.padding = unit(0.3, "lines"), colour= "black", 
                     size= 4, show.legend = F )+ 
    
    
    theme_bw()+
    scale_size(range = c(2.5,6))+
    theme(plot.title = element_text(hjust = 0.5),
          axis.title.x= element_text( size=10, margin = NULL,face="bold"),
          axis.title.y= element_text( size=10, margin = NULL,face="bold"),
          axis.text.x = element_text(face="bold",  size=10),
          axis.text.y = element_text(face="bold",  size=10),
          legend.title = element_blank()) 
  
  return(figure)
}


complex_volcanoF <- function(data, sample = NULL) {
  require(data.table)
  samp <- data[ SAMPLE_ID == sample] %>% as.data.table() #subset for patient
  samp <- samp[!is.na(COMPLEX_LOG2FC)]  #remove missing
  samp <- samp[ !duplicated(samp[, c("SAMPLE_ID", "COMPLEX")]) ]
  
  
  figure <- ggplot(samp, aes(COMPLEX_LOG2FC, -log10(COMPLEX_PVALUE) )) +
    geom_vline(xintercept = 0, color = "grey50") +
    
    geom_point(aes(color= COMPLEX_Z_outlier)) + 
    
    xlab("Protein complex log2(Fold change)") + 
    ylab("-log10(P-value)") +
    scale_color_manual( values = c("gray80", "darkorange"))+
    
    geom_text_repel(data= samp[COMPLEX_Z_outlier == T ], 
                    mapping=aes(COMPLEX_LOG2FC, -log10(COMPLEX_PVALUE), label= COMPLEX), 
                    size= 3,  show.legend = F, min.segment.length = Inf )+ 
    
    theme_bw()+
    theme(plot.title = element_text(hjust = 0.5),
          axis.title.x= element_text( size=10, margin = NULL,face="bold"),
          axis.title.y= element_text( size=10, margin = NULL,face="bold"),
          axis.text.x = element_text(face="bold",  size=10),
          axis.text.y = element_text(face="bold",  size=10),
          legend.title = element_blank(),
          legend.position = "none") 
  
  return(figure)
}


###########################################################
###########################################################

rank_protein <- function(data, gene = NULL) {
  require(data.table)
  samp <- data[ geneID == gene] %>% as.data.table() #subset for patient
  samp <- samp[!is.na(PROTEIN_INT)]  #remove missing

  samp[, outlier_class := "non_outlier"]
  samp[PROTEIN_outlier == T, outlier_class := "protein-only"]
  samp[ is.na(Semantic_sim), Semantic_sim:= 0]
  
  figure <- ggplot(samp, aes(rank_protein, PROTEIN_INT )) +
    
    geom_point(aes(color= outlier_class, shape= gene_class, size = Semantic_sim )) + 
    ggtitle(gene) +
    xlab("Sample rank") + 
    ylab("Normalized protein intensity") +
    scale_color_manual(breaks = c("non_outlier","RNA-only","protein-only", "RNA-and-protein"),
                       values = outlier_colors)+
    
    scale_shape_manual(values = gene_shapes) +
    
    
    
    geom_text_repel(data= samp[outlier_class != "non_outlier" & (gene_class == "no rare variant" | gene_class == "no data") ], 
                    mapping=aes(rank_protein, PROTEIN_INT, label= SAMPLE_ID), colour= "grey20", 
                    size= 2.8, show.legend = F)+ 
    
    geom_text_repel(data= samp[outlier_class != "non_outlier" & gene_class == "1 rare variant" ], 
                    mapping=aes(rank_protein, PROTEIN_INT, label= SAMPLE_ID), colour= "black", 
                    size= 3.3, show.legend = F)+ 
    
    geom_label_repel(data= samp[outlier_class != "non_outlier" & gene_class == "rare pot. biallelic variants"], 
                     mapping=aes(rank_protein, PROTEIN_INT, label= SAMPLE_ID ), 
                     box.padding = unit(0.3, "lines"), colour= "black", 
                     size= 4, show.legend = F )+ 
    
    
    theme_bw()+
    scale_size(range = c(2.5,6))+
    theme(plot.title = element_text(hjust = 0.5, size=14),
          axis.title.x= element_text( size=10, margin = NULL,face="bold"),
          axis.title.y= element_text( size=10, margin = NULL,face="bold"),
          axis.text.x = element_text(face="bold",  size=10),
          axis.text.y = element_text(face="bold",  size=10),
          legend.title = element_blank(),
          plot.margin = margin(0, 0, 0, 0, "cm") ) 
  
  return(figure)
}

rank_rna <- function(data, gene = NULL) {
  require(data.table)
  samp <- data[ geneID == gene] %>% as.data.table() #subset for patient
  samp <- samp[!is.na(RNA_LOG2FC)]  #remove missing
  
  samp[, outlier_class := "non_outlier"]
  samp[RNA_outlier == T, outlier_class := "RNA-only"]
  samp[ is.na(Semantic_sim), Semantic_sim:= 0]
  
  figure <- ggplot(samp, aes(rank_rna, normcounts +1 )) +
    
    geom_point(aes(color= outlier_class, shape= gene_class, size = Semantic_sim )) + 
    ggtitle(gene) +
    xlab("Sample rank") + 
    ylab("Normalized counts + 1") +
    scale_color_manual(breaks = c("non_outlier","RNA-only","protein-only", "RNA-and-protein"),
                       values = outlier_colors)+
    
    scale_shape_manual(values = gene_shapes) +
    
    
    
    geom_text_repel(data= samp[outlier_class != "non_outlier" & (gene_class == "no rare variant" | gene_class == "no data")  ], 
                    mapping=aes(rank_rna, normcounts + 1, label= SAMPLE_ID), colour= "grey20", 
                    size= 2.8, show.legend = F)+ 
    
    geom_text_repel(data= samp[outlier_class != "non_outlier" & gene_class == "1 rare variant" ], 
                    mapping=aes(rank_rna, normcounts +1 , label= SAMPLE_ID), colour= "black", 
                    size= 3.3, show.legend = F)+ 
    
    geom_label_repel(data= samp[outlier_class != "non_outlier" & gene_class == "rare pot. biallelic variants"], 
                     mapping=aes(rank_rna, normcounts +1 , label= SAMPLE_ID ), 
                     box.padding = unit(0.3, "lines"), colour= "black", 
                     size= 4, show.legend = F )+ 
    
    
    theme_bw()+
    scale_size(range = c(2.5,6))+
    theme(plot.title = element_text(hjust = 0.5, size=14),
          axis.title.x= element_text( size=10, margin = NULL,face="bold"),
          axis.title.y= element_text( size=10, margin = NULL,face="bold"),
          axis.text.x = element_text(face="bold",  size=10),
          axis.text.y = element_text(face="bold",  size=10),
          legend.title = element_blank(),
          plot.margin = margin(0, 0, 0, 0, "cm") ) 
  
  return(figure)
}

rank_protein_rna <- function(data, gene = NULL) {
  require(data.table)
  samp <- data[ geneID == gene] %>% as.data.table() #subset for patient
  samp <- samp[!is.na(prot_rna)]  #remove missing
  samp[ is.na(Semantic_sim), Semantic_sim:= 0]
  
  figure <- ggplot(samp, aes(rank_prot_rna, prot_rna )) +
    
    geom_point(aes(color= outlier_class, shape= gene_class, size = Semantic_sim )) + 
    ggtitle(gene) +
    xlab("Sample rank") + 
    ylab("Protein / RNA ratio") +
    scale_color_manual(breaks = c("non_outlier","RNA-only","protein-only", "RNA-and-protein"),
                       values = outlier_colors)+
    
    scale_shape_manual(values = gene_shapes) +
    
    
    
    geom_text_repel(data= samp[outlier_class != "non_outlier" & (gene_class == "no rare variant" | gene_class == "no data") ], 
                    mapping=aes(rank_prot_rna, prot_rna, label= SAMPLE_ID), colour= "grey20", 
                    size= 2.8, show.legend = F)+ 
    
    geom_text_repel(data= samp[outlier_class != "non_outlier" & gene_class == "1 rare variant" ], 
                    mapping=aes(rank_prot_rna, prot_rna, label= SAMPLE_ID), colour= "black", 
                    size= 3.3, show.legend = F)+ 
    
    geom_label_repel(data= samp[outlier_class != "non_outlier" & gene_class == "rare pot. biallelic variants"], 
                     mapping=aes(rank_prot_rna, prot_rna , label= SAMPLE_ID ), 
                     box.padding = unit(0.3, "lines"), colour= "black", 
                     size= 4, show.legend = F )+ 
    
    
    theme_bw()+
    scale_size(range = c(2.5,6))+
    theme(plot.title = element_text(hjust = 0.5, size=14),
          axis.title.x= element_text( size=10, margin = NULL,face="bold"),
          axis.title.y= element_text( size=10, margin = NULL,face="bold"),
          axis.text.x = element_text(face="bold",  size=10),
          axis.text.y = element_text(face="bold",  size=10),
          legend.title = element_blank(),
          plot.margin = margin(0, 0, 0, 0, "cm") ) 
  
  return(figure)
}


###########################################################
###########################################################

# create a theme for plotting
mytheme <- gridExtra::ttheme_default(
  core = list(padding=unit(c(0, 0), "mm"))
)

