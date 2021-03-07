#'---
#' title: Proportion of outliers with variants and phenotype simmilarity
#' author: smirnovd
#' wb:
#'  input:
#'  - config: 'src/config.R'
#'  - patient_omics: '`sm config["PROC_DATA"] + "/integration/patient_omics_full.RDS"`'
#'  output:
#'  - proportions_rare_pb: '`sm config["PROC_DATA"] + "/variant_tables/proportions_rare_biallelic.tsv"`'
#' output: 
#'   html_document:
#'    code_folding: hide
#'    code_download: TRUE
#'---

# Load plotting config and functions
source(snakemake@input$config)
source("src/functions/Integration/integrate_annotate_omics.R")
source("src/functions/variant_enrichment.R")


# Read integrated omics file 
rp <- readRDS(snakemake@input$patient_omics) %>% as.data.table()

# Subset cases with WES and RNA-seq data avaliable 
rp <- rp[WES_avaliable == T & RNA_seq_avaliable == T]
paste("Number of samples with WES export and RNA-seq avaliable:", uniqueN(rp$SAMPLE_ID))


#Keep only genes, with both RNA and protein measured 
rp <- rp[ gene_detected == "RNA and protein detected"]

# Filter for the genes, detected as outliers at least once
outliers <- unique(rp[outlier_class != "non_outlier"]$geneID)
rp <- rp[geneID %in% outliers]


# Add up- down- outlier class
rp <- add_up_down_class(rp, 
                        Padj_threshold =  PADJ_THRESHOLD,
                        Zscore_threshold = ZSCORE_THRESHOLD, 
                        l2FC_threshold = LOG2FC_THRESHOLD)


# Add phenotype semantic simillarity class
rp[, ss_cat:= Semantic_sim > 2]
rp[is.na(ss_cat), ss_cat := F]


# Add variant type annotation
rp[ , var_type := "no rare variant"] 
rp[ rare == T , var_type := "rare"] 
rp[ potential_biallelic == T , var_type := "potential_biallelic"] 
rp[potential_biallelic == T & ss_cat == T,  var_type:= "potential biallelic SemSim > 2"]

# Subset necessary columns
rp <- rp[ , c("SAMPLE_ID", "geneID", "var_type", "up_down_outlier") ]


# Calculate proportions of outliers
pr <- rp[, .N, by = .(up_down_outlier, var_type)]
pr[, total := sum(N), by = .(up_down_outlier)]
pr[, prop := N/sum(N), by = up_down_outlier]
pr[, type := "pb"]


write_tsv(pr, snakemake@output$proportions_rare_pb)

pr[, var_type := factor(var_type, levels = c("no rare variant", "rare", "potential_biallelic", "potential biallelic SemSim > 2"))]
pr[, outlier_class_label := paste0(up_down_outlier ,'\n', "(n = ", total, ")") ]


#' ### Underexpression outliers 
#' Performed on the not representative subset of the data, to illustrate the analysis 
pr <- pr[up_down_outlier %in% c("RNA_underexpression", "Protein_underexpression", "RNA_Protein_underexpression", "non_outlier")]

#+ fig.width=13, fig.height=3
ggplot(pr, aes(outlier_class_label, prop)) +
  geom_bar(stat= 'identity', aes(fill = var_type)) +
  scale_fill_manual(values = variant_colors2 ) +
  coord_flip(ylim = c(0,1)) +
  scale_y_continuous(breaks=seq(0,1, 0.1), labels=scales::percent) +
  labs( y = "Proportion of outliers with rare variants")+
  theme_classic()+
  theme(legend.position="top",  axis.title.y = element_blank() ,legend.direction = "horizontal",
        legend.title = element_blank(),
        axis.text.y = element_text(face="bold", size=12)) +
  guides(fill = guide_legend(nrow = 1))

