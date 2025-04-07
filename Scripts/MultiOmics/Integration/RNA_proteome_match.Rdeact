#'---
#' title: RNAseq - proteome matching and Fig S4
#' author: vyepez
#' wb:
#'  input: 
#'  - sample_annotation: '`sm config["ANNOTATION"]`'
#'  - protein_gene_mat: '`sm config["RAW_Protein"]`'
#'  - ods_fib: '`sm config["PROC_DATA"] + "/outrider/ods.Rds"`'
#'  - rna_prot_cor: '`sm config["PROC_DATA"] + "/integration/rna_protein_cor.tsv"`'
#'  output:
#'  - fig1: '`sm config["FIGURE_DIR"] + "/Supplementary_figures/S_Fig4_a.pdf"`'
#'  - fig2: '`sm config["FIGURE_DIR"] + "/Supplementary_figures/S_Fig4_b.pdf"`'
#' output: 
#'   html_document:
#'    code_folding: hide
#'    code_download: TRUE
#'---


suppressPackageStartupMessages({
  library(data.table)
  library(magrittr)
  library(ggplot2)
  library(ggthemes)
  library(dplyr)
  library(matrixStats)
  library(OUTRIDER)
  library(mclust)
  library(LSD)
  library(plotly)
  library(pheatmap)
  library(ggplotify)
})



# Read sample annotation and subset to samples that we have either RNA or Proteome
sa <- fread(snakemake@input$sample_annotation)

# Prepare annotation 
sa <- sa[ !(NORMALIZATION_SAMPLE == T & USE_FOR_PROTEOMICS_PAPER == F)]
sa[, rna := tstrsplit(SAMPLE_ID, ".2", fixed=TRUE)]
sa[TREATMENT == "HEAT", rna := tstrsplit(SAMPLE_ID, "_heat", fixed=TRUE)]
sa[TREATMENT == "GAL", rna := tstrsplit(SAMPLE_ID, "_GAL", fixed=TRUE)]
sa[TREATMENT == "DNAJC30", rna := tstrsplit(SAMPLE_ID, "_T_DNAJC30", fixed=TRUE)]
sa[TREATMENT == "MRPL38", rna := tstrsplit(SAMPLE_ID, "_T_MRPL38", fixed=TRUE)]
sa <- sa[is.na(TREATMENT) & TISSUE == "FIBROBLAST"]
sa <- sa[USE_FOR_PROTEOMICS_PAPER == T]
sa[USE_FOR_PROTEOMICS_PAPER == T, rna :=  SAMPLE_ID]

sa[ , RNA_ID := paste0(rna, "R")]
sa[ , PROTEOME_ID := paste0(SAMPLE_ID, "P")]
sa <- sa[! (is.na(RNA_ID) & is.na(PROTEOME_ID))]

# Read protein matrix, subset, log-transform and center
protein_gene_mat <- fread(snakemake@input$protein_gene_mat) %>% as.data.frame()
rownames(protein_gene_mat) <- protein_gene_mat$geneID
protein_gene_mat$geneID <- NULL
setnames(protein_gene_mat, sa$SAMPLE_ID, sa$PROTEOME_ID, skip_absent=TRUE)
protein_gene_mat <- as.matrix(protein_gene_mat)
protein_gene_mat[protein_gene_mat < 1e5] <- 0


rna_prot_dt <- sa[PROTEOME_ID %in% colnames(protein_gene_mat)  & TISSUE == 'FIBROBLAST', .(PROTEOME_ID, RNA_ID)] 


# Subset Proteomes
protein_gene_mat <- protein_gene_mat[, rna_prot_dt$PROTEOME_ID]
protein_gene_mat[protein_gene_mat < 1e5] <- NA
protein_gene_mat <- log(protein_gene_mat + 1)
protein_gene_mat <- protein_gene_mat - rowMeans2(protein_gene_mat, na.rm = T)
row.names(protein_gene_mat) <- toupper(row.names(protein_gene_mat))


# Read RNA data
ods <- readRDS(snakemake@input$ods)
counts <- counts(ods, normalized = F)
counts <- t(t(counts) / sizeFactors(ods))
counts[counts < 50] <- NA
counts <- log(counts + 1) - rowMeans2(log(counts + 1), na.rm = T)

counts <- as.data.frame(counts)
setnames(counts, sa$SAMPLE_ID, sa$RNA_ID, skip_absent=TRUE)
counts <- as.matrix(counts)

rna_prot_dt <- rna_prot_dt[RNA_ID %in% colnames(counts) ]
rna_prot_dt <- rna_prot_dt[order(RNA_ID, decreasing = T)]

#' Dimensions of RNA count table
counts <- counts[, unique(rna_prot_dt$RNA_ID)]
dim(counts)

#' Dimensions of protein intensity table
protein_gene_mat <- protein_gene_mat[, rna_prot_dt$PROTEOME_ID]
dim(protein_gene_mat)


# Read rna - protein correlation data 
rpc <- fread(snakemake@input$rna_prot_cor)
# Genes with significantly correlated RNA and protein levels
rpc_sign <- rpc[p.value < 0.05  ] # significant == T

common_genes <- intersect(row.names(counts), row.names(protein_gene_mat)) 
disp <- dispersions(ods)
names(disp) <- rownames(ods)
disp <- disp[ common_genes]
top_disp_genes <- names(head(sort(disp, decreasing = T), 1000))
common_genes <- intersect(common_genes, top_disp_genes)

# Uncomment for the additional filtration for the genes with signif correlated RNA and protein levels
# common_genes <- intersect(common_genes, rpc_sign$geneID)

#' Number of common genes
length(common_genes)

#' Dimensions of RNA count table after the gene subset
cm <- counts[common_genes, ]
dim(cm)

#' Dimensions of protein intensity table  after the gene subset
pm <- protein_gene_mat[common_genes, ]
dim(pm)



#' Correlation matrix 
x = sapply(1:ncol(cm), function(i){
  sapply(1:ncol(pm),
         function(j) cor.test(pm[,j], cm[,i], method = 'spearman')$estimate)
}
)
hist(x, main = 'Correlation of all RNA - Protein Permutations'); abline(v = 0, col = 'red', lty = 'dashed')


colnames(x) <- colnames(cm)
rownames(x) <- colnames(pm)
#apply(x, 1, which.max)

#' Number of times where the highest correlation was not on the annotated RNA
sum(apply(x, 1, which.max) != 1:nrow(x))

#' # Correlation of all RNA - Protein combinations
#+ fig.width=18, fig.height=18
Fig_S4a <- as.ggplot(pheatmap(x, cluster_rows=F, cluster_cols=F, 
                              main = "Transcriptome - proteome correlation per sample"))



#' Save supplementary figure 4a
pdf(snakemake@output$fig1,  
    width = 15, height = 15,  useDingbats=FALSE )
print(Fig_S4a) 
dev.off()


# y <- sapply(1:ncol(pm), function(j) cor.test(pm[,j], cm[,j], method = 'spearman')$estimate)
# names(y) <- colnames(pm)


y <- sapply(1:nrow(rna_prot_dt), function(j) cor.test(pm[,rna_prot_dt[j,]$PROTEOME_ID ], cm[,rna_prot_dt[j,]$RNA_ID], method = 'spearman')$estimate)
names(y) <- rna_prot_dt$PROTEOME_ID 

#' ## Expectation Maximation to see the classes' separation
mod <- Mclust(y)
mod4 <- densityMclust(y)
#' Number of groups the algorithm suggests
mod4$G

#' Proteome ids that belong to other cluster
names(mod$classification[mod$classification == 1]) %>% sort

#' Compute the best separation between the classes
mus <- mod4$parameters$mean
sigmas <- sqrt(mod4$parameters$variance$sigmasq)
fr <- function(x) {
  1/sqrt(2*pi*sigmas[1]^2) * exp(-(x-mus[1])^2/(2*sigmas[1]^2)) + 1/sqrt(2*pi*sigmas[2]^2) * exp(-(x-mus[2])^2/(2*sigmas[2]^2))
}
op <- optim(.3, fr, method = 'Brent', lower = -.5, upper = .8)

plot(mod4, what = "density", data = y, breaks = 15, xlab = "Correlation of RNA - Protein annotated samples")
abline(v = op$par, col = 'red', lty = 'dashed')


#' ## Find possible matches
dt <- as.data.table(melt(x))
setnames(dt, old = c("Var1", "Var2"), c("Proteome_ID", "RNA_ID"))
dt[, max_corr := value == max(value), by = Proteome_ID]
dt[, aux := paste(Proteome_ID, RNA_ID, sep = "-")]
rna_prot_dt[, aux := paste(PROTEOME_ID, RNA_ID, sep = "-")]
dt[, right_annot := aux %in% rna_prot_dt$aux]
dt[, aux := NULL]
DT::datatable(dt[Proteome_ID %in% names(mod$classification[mod$classification == 1]), ][max_corr == T | right_annot == T],  options = list(pageLength = 20))



# Create separate dataset wo confirmed matches
dtX <- dt[max_corr == T & right_annot == T]
dts <- dt[! ( (Proteome_ID %in% dtX$Proteome_ID) | (RNA_ID %in% dtX$RNA_ID) ) ]
dts[, max_corr := value == max(value), by = Proteome_ID]
# View(dts[right_annot == T |max_corr == T])



#' ## Plot all correlations of proteomes 
dt[, Proteome_ID := as.character(Proteome_ID)]
plotlist = list()
for(s in unique(dt$Proteome_ID)){
  g <- ggplot(dt[Proteome_ID == s], aes(reorder(RNA_ID, value), value)) + geom_point(aes(col = right_annot)) + theme_bw() + 
    theme(axis.text.x=element_blank()) + scale_color_calc() + 
    labs(title = s, x = 'Rank', y = 'RNA - Protein Correlation')
  plotlist[[s]] = ggplotly(g)
}

#' Uncomment to produce the full list of figures
#' May cause accidental browser crash 
# htmltools::tagList(plotlist)

#' Supplementary figure 4b
examp <- ggplot(dt[Proteome_ID == "OM35261P"], aes(reorder(RNA_ID, value), value)) + geom_point(aes(col = right_annot)) + theme_bw() + 
  theme(axis.text.x=element_blank()) + scale_color_calc() + 
  labs(title = "OM35261P", x = 'Rank', y = 'RNA - Protein Correlation')

#+ fig.width=7, fig.height=5
ggplotly(examp)


pdf(snakemake@output$fig2,  
    width = 7, height = 5,  useDingbats=FALSE )
plot(examp)
dev.off()

