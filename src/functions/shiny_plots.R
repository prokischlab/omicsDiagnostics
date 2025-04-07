library(plotly)



plot_patient_integration_plotly <- function(data) {
  p <- ggplot(data[ !is.na(PROTEIN_ZSCORE ) & !is.na(RNA_ZSCORE )  ], aes(RNA_ZSCORE, PROTEIN_ZSCORE)) +
    geom_vline(xintercept = 0, color = "grey50") +
    geom_hline(yintercept = 0, color = "grey50") +
    geom_abline(color = "grey60") +
    geom_point(aes(color = outlier_class, shape = gene_class, size = SemSim, key = sample_gene  )) +
    xlab("RNA Z-score") +
    ylab("Protein Z-score") +
    scale_color_manual(breaks = c("non_outlier", "RNA_outlier", "Protein_outlier", "RNA_Protein_outlier"), values = outlier_colors) +
    scale_shape_manual(values = gene_shapes) +
    theme_bw() +
    theme(axis.title.x = element_text(size = 9, margin = NULL, face = "bold"),
          axis.title.y = element_text(size = 9, margin = NULL, face = "bold"),
          legend.position = "bottom")
  
  ggplotly(p, tooltip = c( "sample_gene", "key", "y", "x", "size", "shape" )) %>%
    layout(legend = list(orientation = "h", y = -0.2))
}

plot_protein_volcano_plotly <- function(data) {
  p <- ggplot(data, aes(PROTEIN_LOG2FC, -log10(PROTEIN_PVALUE))) +
    geom_vline(xintercept = 0, color = "grey50") +
    geom_point(aes(color = outlier_class, shape = gene_class, size = SemSim, key = sample_gene )) +
    xlab("Protein log2(Fold change)") +
    ylab("-log10(P-value)") +
    scale_color_manual(breaks = c("non_outlier", "RNA-only", "protein-only", "RNA-and-protein"), values = outlier_colors) +
    scale_shape_manual(values = gene_shapes) +
    theme_bw() +
    theme(axis.title.x = element_text(size = 9, margin = NULL, face = "bold"),
          axis.title.y = element_text(size = 9, margin = NULL, face = "bold"),
          legend.position = "bottom")

  
  ggplotly(p, tooltip = c( "sample_gene", "key", "y", "x", "size", "shape" )) %>%
    layout(legend = list(orientation = "h", y = -0.2))
}



plot_rna_volcano_plotly <- function(data) {
  p <- ggplot(data, aes(RNA_LOG2FC, -log10(RNA_PVALUE))) +
    geom_vline(xintercept = 0, color = "grey50") +
    geom_point(aes(color = outlier_class, shape = gene_class, size = SemSim, key = sample_gene)) +
    xlab("RNA log2(Fold change)") +
    ylab("-log10(P-value)") +
    scale_color_manual(breaks = c("non_outlier", "RNA-only", "protein-only", "RNA-and-protein"), values = outlier_colors) +
    scale_shape_manual(values = gene_shapes) +
    theme_bw() +
    theme(axis.title.x = element_text(size = 9, margin = NULL, face = "bold"),
          axis.title.y = element_text(size = 9, margin = NULL, face = "bold"),
          legend.position = "bottom")
  
  ggplotly(p, tooltip = c("sample_gene", "key", "y", "x", "size", "shape")) %>%
    layout(legend = list(orientation = "h", y = -0.2))
}

plot_anatogram <- function(sample_id , patients_data, sa_data) {
  if (!is.na(sample_id)) {
    dat <- patients_data[SAMPLE_ID == sample_id]
  } else {
    dat <- patients_data
  }
  
  
  plot_org <- gganatogram(
    data = patients_data[SAMPLE_ID == sample_id], 
    fillOutline = 'white', 
    organism = 'human', 
    sex = unique(sa_data[SAMPLE_ID == sample_id]$gender), 
    fill = "colour"
  ) + 
    theme_void() 
  return(plot_org)
}


plot_complex_volcano_plotly <- function(sample_id , complex_data) {
  data <- complex_data[SAMPLE_ID == sample_id]
  p <- ggplot(data, aes(COMPLEX_LOG2FC, -log10(COMPLEX_PVALUE))) +
    geom_vline(xintercept = 0, color = "grey50") +
    geom_point(aes(color = COMPLEX_Z_outlier, key = COMPLEX)) +
    xlab("Protein complex log2(Fold change)") +
    ylab("-log10(P-value)") +
    ggtitle(paste0("Aberrant protein complexes for ", sample_id  ) ) +
    scale_color_manual(values = c("gray80", "darkorange")) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5),
          axis.title.x = element_text(size = 10, margin = NULL, face = "bold"),
          axis.title.y = element_text(size = 10, margin = NULL, face = "bold"),
          axis.text.x = element_text(face = "bold", size = 10),
          axis.text.y = element_text(face = "bold", size = 10),
          legend.title = element_blank(),
          legend.position = "none")
  
  ggplotly(p, tooltip = c("key", "COMPLEX", "COMPLEX_LOG2FC", "-log10(COMPLEX_PVALUE)")) %>%
    layout(legend = list(orientation = "h", y = -0.2))
}


plot_rank_rna_plotly <- function(data, gene) {
  require(data.table)
  samp <- data[geneID == gene] %>% as.data.table()  # subset for patient
  samp <- samp[!is.na(normcounts)]  # remove missing
  
  samp[, outlier_class := "non_outlier"]
  samp[RNA_outlier == T, outlier_class := "RNA-only"]
  samp[is.na(SemSim), SemSim := 0]
  
  p <- ggplot(samp, aes(rank_rna, normcounts +1)) +
    geom_point(aes(color = outlier_class, shape = gene_class, size = SemSim, key = SAMPLE_ID)) +
    ggtitle(paste0("Rank of ", gene , " gene expression") ) +
    xlab("Sample rank") +
    ylab("Normalized protein intensity") +
    scale_color_manual(breaks = c("non_outlier", "RNA-only", "protein-only", "RNA-and-protein"), values = outlier_colors) +
    scale_shape_manual(values = gene_shapes) +
    theme_bw() +
    theme(axis.title.x = element_text(size = 10, margin = NULL, face = "bold"),
          axis.title.y = element_text(size = 10, margin = NULL, face = "bold"),
          axis.text.x = element_text(face = "bold", size = 10),
          axis.text.y = element_text(face = "bold", size = 10),
          legend.title = element_blank(), legend.position = "none")
  
  ggplotly(p, tooltip = c("key", "rank_rna", "normcounts", "size", "color", "shape")) # %>% layout(legend = list(orientation = "h", y = -0.2))
}



plot_rank_protein_plotly <- function(data, gene) {
  require(data.table)
  samp <- data[geneID == gene] %>% as.data.table()  # subset for patient
  samp <- samp[!is.na(PROTEIN_INT)]  # remove missing
  
  samp[, outlier_class := "non_outlier"]
  samp[PROTEIN_outlier == T, outlier_class := "protein-only"]
  samp[is.na(SemSim), SemSim := 0]
  
  p <- ggplot(samp, aes(rank_protein, PROTEIN_INT)) +
    geom_point(aes(color = outlier_class, shape = gene_class, size = SemSim, key = SAMPLE_ID)) +
    ggtitle(paste0("Rank of ", gene , " protein intensity") ) +
    xlab("Sample rank") +
    ylab("Normalized protein intensity") +
    scale_color_manual(breaks = c("non_outlier", "RNA-only", "protein-only", "RNA-and-protein"), values = outlier_colors) +
    scale_shape_manual(values = gene_shapes) +
    theme_bw() +
    theme(axis.title.x = element_text(size = 10, margin = NULL, face = "bold"),
          axis.title.y = element_text(size = 10, margin = NULL, face = "bold"),
          axis.text.x = element_text(face = "bold", size = 10),
          axis.text.y = element_text(face = "bold", size = 10),
          legend.title = element_blank(), legend.position = "none")
  
  ggplotly(p, tooltip = c("key", "rank_protein", "PROTEIN_INT", "size", "color", "shape")) # %>% layout(legend = list(orientation = "h", y = -0.2))
}
