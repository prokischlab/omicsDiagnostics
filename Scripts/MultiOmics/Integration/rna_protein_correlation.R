#'---
#' title: RNA vs protein correlation  
#' author: smirnovd
#' wb:
#'  input:
#'  - sample_annotation: '`sm config["ANNOTATION"]`'
#'  - disease_genes: '`sm config["DATASETS"] + "/disease_genes.tsv"`'
#'  - patient_omics: '`sm config["PROC_DATA"] + "/integration/patient_omics.RDS"`'
#'  output:
#'  - rna_prot_cor: '`sm config["PROC_DATA"] + "/integration/rna_protein_cor.tsv"`'
#' output: 
#'   html_document:
#'    code_folding: hide
#'    code_download: TRUE
#'---

# Load plotting config and functions
source('src/config.R')

# READ ANNOTATION
# sa <- fread('/s/project/mitoMultiOmics/multiOMICs_integration/raw_data/proteomics_annotation.tsv')
sa <- fread(snakemake@input$sample_annotation)
sa <- sa[USE_FOR_PROTEOMICS_PAPER == T]



# Load disease genes table
# dis_genes <- fread('/s/project/mitoMultiOmics/multiOMICs_integration/datasets/disease_genes.tsv')
dis_genes <- fread(snakemake@input$disease_genes)


# Read integrated omics file 
# rp <- readRDS("/s/project/mitoMultiOmics/multiOMICs_integration/processed_data/integration/patient_omics.RDS") %>% as.data.table()
rp <- readRDS(snakemake@input$patient_omics) %>% as.data.table()
rp <- rp[ gene_detected == "RNA and protein detected"]


##################################################

#' # Z-scores
range <- c(rp$RNA_ZSCORE, rp$PROTEIN_ZSCORE)
#+ fig.width=7, fig.height=7
p1 <- ggplot(rp, aes(RNA_ZSCORE, PROTEIN_ZSCORE) )+
  geom_hex( )+
  stat_cor( aes(label = paste(..r.label.. , ..rr.label.., ..p.label.., sep = "~`,`~")), method = "pearson", na.rm = T)+
  geom_smooth(method = lm, se = FALSE)+
  coord_fixed( xlim = c( min(range , na.rm = T) , max( range  , na.rm = T) ),  
               ylim = c( min(range , na.rm = T),  max( range  , na.rm = T)  )) + 
  xlab("RNA zScore") + 
  ylab("Protein zScore") +
  theme_bw()+
  theme(plot.title = element_text(size=14, hjust = 0.5, face="bold"),
        axis.title.x= element_text( size=10, margin = NULL,face="bold"),
        axis.title.y= element_text( size=10, margin = NULL,face="bold"),
        axis.text.x= element_text( size=9, margin = NULL,face="bold"),
        axis.text.y= element_text( size=9, margin = NULL,face="bold"))
p1
##################################################

#' # log2(Fold change)
rp$RNA_LOG2FC <- log2(rp$RNA_FC)
rp$PROTEIN_LOG2FC <- log2(rp$PROTEIN_FC)

rpx <- rp[!is.na(RNA_LOG2FC ) & !is.na(PROTEIN_LOG2FC )]
rpx <- rp[!is.infinite(RNA_LOG2FC ) & !is.infinite(PROTEIN_LOG2FC )]

range <- c(rpx$RNA_LOG2FC , rpx$PROTEIN_LOG2FC)
#+ fig.width=7, fig.height=7
ggplot(rpx, aes(RNA_LOG2FC, PROTEIN_LOG2FC ) )+
  geom_hex( )+
  geom_vline(xintercept = 0, color = "grey50") +
  geom_hline(yintercept = 0, color = "grey50") +
  stat_cor( aes(label = paste(..r.label.. , ..rr.label.., ..p.label.., sep = "~`,`~")), method = "pearson", na.rm = T)+
  geom_smooth(method = lm, se = FALSE)+
  coord_fixed( xlim = c( min(range , na.rm = T) , max( range  , na.rm = T) ),  
               ylim = c( min(range , na.rm = T),  max( range  , na.rm = T)  )) + 
  xlab("RNA log2(fold change)") + 
  ylab("Protein log2(fold change)") +
  theme_bw()+
  theme(plot.title = element_text(size=14, hjust = 0.5, face="bold"),
        axis.title.x= element_text( size=10, margin = NULL,face="bold"),
        axis.title.y= element_text( size=10, margin = NULL,face="bold"),
        axis.text.x= element_text( size=9, margin = NULL,face="bold"),
        axis.text.y= element_text( size=9, margin = NULL,face="bold"))

##################################################
#' # Abundance
rp$log10_RNA_counts <- log10(rp$normcounts)
rp$log10_Protein_intensity <- log10(rp$PROTEIN_INT)

rpx <- rp[!is.na(log10_RNA_counts ) & !is.na(log10_Protein_intensity )]
rpx <- rp[!is.infinite(log10_RNA_counts ) & !is.infinite(log10_Protein_intensity )]

range <- c(rpx$log10_RNA_counts , rpx$log10_Protein_intensity)
#+ fig.width=7, fig.height=7
p2 <- ggplot(rpx, aes(log10_RNA_counts, log10_Protein_intensity ) )+
  geom_hex( )+
  geom_vline(xintercept = 0, color = "grey50") +
  geom_hline(yintercept = 0, color = "grey50") +
  stat_cor( aes(label = paste(..r.label.. , ..rr.label.., ..p.label.., sep = "~`,`~")), method = "pearson", na.rm = T)+
  geom_smooth(method = lm, se = FALSE)+
  coord_fixed( xlim = c( min(range , na.rm = T) , max( range  , na.rm = T) ),  
               ylim = c( min(range , na.rm = T),  max( range  , na.rm = T)  )) + 
  xlab("RNA log10(counts)") + 
  ylab("Protein log10(intensity)") +
  theme_bw()+
  theme(plot.title = element_text(size=14, hjust = 0.5, face="bold"),
        axis.title.x= element_text( size=10, margin = NULL,face="bold"),
        axis.title.y= element_text( size=10, margin = NULL,face="bold"),
        axis.text.x= element_text( size=9, margin = NULL,face="bold"),
        axis.text.y= element_text( size=9, margin = NULL,face="bold"))

p2

##################################################

#' # Gene - protein correlation Spearman

gp_s <- as.data.table(rpx[ , cor.test(log10_RNA_counts, log10_Protein_intensity, method="spearman")[-2], by= geneID])
setnames(gp_s, "estimate", "rho" )
med <- median(gp_s$rho)
gp_s[ , padj := p.adjust(p.value, method = "BH")]
gp_s[, significant := padj < 0.05]




# Save results
write_tsv(gp_s, snakemake@output$rna_prot_cor)
# write_tsv(gp_s, "/s/project/mitoMultiOmics/multiOMICs_integration/processed_data/integration/rna_protein_cor.tsv")



#+ fig.width=8, fig.height=6
p3 <- ggplot(gp_s, aes(rho, fill = significant))+
  geom_histogram(bins = 90)+
  theme_classic()+
  geom_vline(xintercept = 0, color = "grey50") +
  geom_vline(xintercept = med, linetype = "dashed") +
  scale_fill_brewer()+
  coord_cartesian( xlim = c( -1 , 1 )) + 
  ylab("# of genes") + 
  xlab("Spearman's rho") + 
  ggtitle("Gene - protein correlation") +
  theme(plot.title = element_text(size=14, hjust = 0.5, face="bold"),
        axis.title.x= element_text( size=12, margin = NULL,face="bold"),
        axis.title.y= element_text( size=12, margin = NULL,face="bold"),
        axis.text.x= element_text( size=10, margin = NULL,face="bold"),
        axis.text.y= element_text( size=10, margin = NULL,face="bold"))

#' ## Disease genes 

gp_s_dg <- merge(dis_genes , gp_s, by = "geneID")
#+ fig.width=9, fig.height=7
p4 <- ggplot(gp_s_dg, aes(rho, fill = significant))+
  geom_histogram()+
  theme_classic()+
  geom_vline(xintercept = 0, color = "grey50") +
  geom_vline(xintercept = med, linetype = "dashed") +
  scale_fill_brewer()+
  #coord_cartesian( xlim = c( -1 , 1 )) + 
  ylab("# of genes") + 
  xlab("Spearman's rho") + 
  facet_wrap( ~DISEASE, scales = "free")+
  theme(plot.title = element_text(size=14, hjust = 0.5, face="bold"),
        axis.title.x= element_text( size=12, margin = NULL,face="bold"),
        axis.title.y= element_text( size=12, margin = NULL,face="bold"),
        axis.text.x= element_text( size=10, margin = NULL,face="bold"),
        axis.text.y= element_text( size=10, margin = NULL,face="bold"))


#' ## Causal genes
causal_genes <- gp_s[ geneID %in% sa[CATEGORY %in% c("I", "IIa", "III")]$KNOWN_MUTATION, c("geneID", "rho", "p.value", "padj",  "significant")]

annot <- rp[causal_gene == T , c("geneID", "outlier_class")]

causal_genes <- merge(causal_genes, annot, by = "geneID")


ggplot(causal_genes[ outlier_class != "non_outlier"], aes(rho, -log10(p.value)))+
  geom_point( aes(color = outlier_class ) )+
  geom_vline(xintercept = 0, color = "grey50") +
  geom_hline(yintercept = -log10(0.05), color = "grey50", linetype = "dashed") +
  geom_text_repel(aes(label = geneID), size = 2.5)+
  scale_color_manual(breaks = c("non_outlier","RNA-only", "protein-only", "RNA-and-protein"),
                     values = outlier_colors)+
  xlab("Spearman's rho") + 
  theme_bw()

#' ## Genes detected as outliers

annot <- rp[outlier_class != "non_outlier" , c("geneID", "outlier_class")]

#' consider only genes, that were detected as outliers of one class
annot <- annot[ !( geneID %in%  annot[duplicated(geneID) ]$geneID ) ]
outlier_genes <- merge(gp_s, annot, by = "geneID")

outlier_genes$outlier_class <- factor(outlier_genes$outlier_class , levels =  c("RNA-only","protein-only", "RNA-and-protein"))

ggplot(outlier_genes, aes(outlier_class, rho, fill = outlier_class))+
  geom_boxplot()+
  scale_fill_manual(breaks = c("non_outlier","RNA-only","protein-only", "RNA-and-protein"),
                    values = outlier_colors)+
  theme_bw() + 
  ylab("Spearman's rho") + 
  stat_compare_means(comparisons = list(c("RNA-only" ,"protein-only" ),
                                        c("protein-only",  "RNA-and-protein"),
                                        c("RNA-only", "RNA-and-protein") ) )+
  #stat_compare_means(label.y = 1.2) + 
  ggtitle("Genes detected as outliers of one category")


#+ fig.width=14, fig.height=9
(p2 | p1) / (p3 | p4)
