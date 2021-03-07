#'---
#' title: Enrichment of rare variants in outliers  
#' author: smirnovd
#' wb:
#'  input:
#'  - patient_omics: '`sm config["PROC_DATA"] + "/integration/patient_omics_full.RDS"`'
#'  output:
#'  - enrichment_rare_pb: '`sm config["PROC_DATA"] + "/variant_tables/enrichment_rare_biallelic.tsv"`'
#' output: 
#'   html_document:
#'    code_folding: hide
#'    code_download: TRUE
#'---

# Load plotting config and functions
source('src/config.R')
source("src/functions/Integration/integrate_annotate_omics.R")
source("src/functions/variant_enrichment.R")


# Read integrated omics file 
# rp <- readRDS("/s/project/mitoMultiOmics/multiOMICs_integration/processed_data/integration/patient_omics_full.RDS") %>% as.data.table()
rp <- readRDS(snakemake@input$patient_omics) %>% as.data.table()

# Subset cases with WES and RNA-seq data avaliable 
rp <- rp[WES_avaliable == T & RNA_seq_avaliable == T]


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


###################################################
# Subset necessary columns
# Add variant type annotation
rp[ , var_type := "No rare variant"] 
rp[ rare == T , var_type := "Rare variant"] 
rp[ potential_biallelic == T , var_type := "Potential biallelic\n rare variant"] 


# Subset necessary columns
rp <- rp[ , c("SAMPLE_ID", "geneID", "var_type", "outlier_class", "up_down_outlier") ]
rp <- rp[!duplicated(rp), ]


###################################################

# As soon as it is not allowed to share genetic data
# Shared data doesn't allow to reproduce original figure 
# Code below was used to produse figure on the full datset 


# Remove variant annotation 
# rp <- rp[ , c("SAMPLE_ID", "geneID",  "outlier_class", "up_down_outlier") ]
# rp <- rp[!duplicated(rp), ]
# 
# Read variant annotation 
# vt <- fread("../rare_variants_pb.tsv")
# rp <- rp[SAMPLE_ID %in% unique(vt$SAMPLE_ID) ]
# vt <- vt[SAMPLE_ID %in% unique(rp$SAMPLE_ID) ]
# vt <- vt[ geneID %in% unique(rp$geneID)]
# 
# rp <- merge(vt, rp, by = c("SAMPLE_ID", "geneID"), all.y = T )
# rp[is.na(var_type), var_type := "No rare variant"]



###################################################



paste("Number of samples with WES export and RNA-seq avaliable:", uniqueN(rp$SAMPLE_ID))
rp[ , Ncases_per_gene := .N, by = geneID]
paste("Minimal number of individuals per gene:" , min(rp$Ncases_per_gene))
# Subset necessary columns
rp <- rp[ , c("SAMPLE_ID", "geneID", "var_type", "outlier_class", "up_down_outlier") ]
rp <- rp[!duplicated(rp), ]




####################

#' # Outlier class
os <- rp[, .N, by = outlier_class]

#' ### Number of outliers per class
DT::datatable(os, caption = "Outlier counts",   escape = F)


# All outliers
enr <- reshape2::dcast(data = rp, geneID +  SAMPLE_ID +  var_type ~ outlier_class , 
                       fun.aggregate = cat.agg)

enrichments <- data.frame()
for (outl in unique(rp$outlier_class)){
  tem_df <- enr
  setnames(tem_df, outl, "outliers")
  tem_df <- tem_df[, c("geneID", "SAMPLE_ID", "outliers", "var_type" )]
  feat <- unique(tem_df$var_type)
  tem_df <- dcast(data = tem_df, geneID + SAMPLE_ID  + outliers ~ var_type, fun.aggregate = cat.agg) 
  enrichment <- enrich(tem_df, outlier= "outliers" ,features= feat)
  enrichment$ci_left <- enrichment$Estim - 1.96 * enrichment$Std
  enrichment$ci_right <- enrichment$Estim + 1.96 * enrichment$Std
  enrichment$outlier_class <- rep(outl, nrow(enrichment))
  enrichments <- rbind(enrichments, enrichment)
  enrichment <- NULL
  enr$outliers <- NULL
}

enrichments <- as.data.table(enrichments)
enrichments[, significant := 1 ]
enrichments[enrichments$Pval >= 0.05 , significant := 0.9 ]
nodelist <- levels(enrichments$Cat)
enrichments$Cat <- factor(enrichments$Cat, levels= c("No rare variant", "Rare variant", "Potential biallelic\n rare variant"))


#' ## Outlier class
#+ fig.width=6, fig.height=9
ggplot(data = enrichments, aes(x = Cat, y = Estim)) +
  geom_pointrange(aes(x = Cat, ymin = ci_left, ymax = ci_right, colour = outlier_class), position = position_dodge(width = 1)) +
  geom_hline(yintercept = 0) +
  geom_vline(color = "gray80", linetype = "dashed", xintercept = which(nodelist=="No rare variant")-0.5) +
  geom_vline(color = "gray80", linetype = "dashed", xintercept = which(nodelist=="Rare variant")-0.5) +
  theme_bw() + 
  xlab('') + 
  ylab('Log odds ratio')  + 
  scale_alpha(guide = 'none')+
  scale_color_manual(values = outlier_colors)+ 
  coord_flip(ylim = c(min(enrichments[enrichments$Pval <0.05, 'ci_left' ] , na.rm = T),  max(enrichments[enrichments$Pval <0.05, "ci_right" ], na.rm = T))) 


# write_tsv(enrichments, "/s/project/mitoMultiOmics/multiOMICs_integration/processed_data/variant_tables/enrichment_rare_biallelic.tsv")
write_tsv(enrichments,  snakemake@output$enrichment_rare_pb)


#####################################
#' # Up- down- regulation outliers
os <- rp[, .N, by = up_down_outlier]


#' ### Number of outliers per class
DT::datatable(os, caption = "Outlier counts",   escape = F)

enr <- reshape2::dcast(data = rp, geneID +  SAMPLE_ID +  var_type ~ up_down_outlier , 
                       fun.aggregate = cat.agg)

enrichments <- data.frame()
for (outl in unique(rp$up_down_outlier)){
  tem_df <- enr
  setnames(tem_df, outl, "outliers")
  tem_df <- tem_df[, c("geneID", "SAMPLE_ID", "outliers", "var_type" )]
  feat <- unique(tem_df$var_type)
  tem_df <- dcast(data = tem_df, geneID + SAMPLE_ID  + outliers ~ var_type, fun.aggregate = cat.agg) 
  enrichment <- enrich(tem_df, outlier= "outliers" ,features= feat)
  enrichment$ci_left <- enrichment$Estim - 1.96 * enrichment$Std
  enrichment$ci_right <- enrichment$Estim + 1.96 * enrichment$Std
  enrichment$outlier_class <- rep(outl, nrow(enrichment))
  enrichments <- rbind(enrichments, enrichment)
  enr$outliers <- NULL
}

enrichments <- as.data.table(enrichments)
enrichments[, significant := 1 ]
enrichments[enrichments$Pval >= 0.05 , significant := 0.9 ]
nodelist <- levels(enrichments$Cat)
enrichments$Cat <- factor(enrichments$Cat, levels= c("No rare variant", "Rare variant", "Potential biallelic\n rare variant"))


#' ## Outlier class
#+ fig.width=6, fig.height=9
ggplot(data = enrichments[outlier_class %in% c("RNA_underexpression", "Protein_underexpression", "RNA_Protein_underexpression", "non_outlier")], aes(x = Cat, y = Estim)) +
  geom_pointrange(aes(x = Cat, ymin = ci_left, ymax = ci_right, colour = outlier_class, alpha = significant), position = position_dodge(width = 1)) +
  geom_hline(yintercept = 0) +
  geom_vline(color = "gray80", linetype = "dashed", xintercept = which(nodelist=="No rare variant")-0.5) +
  geom_vline(color = "gray80", linetype = "dashed", xintercept = which(nodelist=="Rare variant")-0.5) +
  theme_bw() + 
  xlab('') + 
  ylab('Log odds ratio')  + 
  scale_alpha(guide = 'none')+
  scale_color_manual(values = outlier_colors)+ 
  coord_flip(ylim = c(min(enrichments[enrichments$Pval <0.05, 'ci_left' ] , na.rm = T),  max(enrichments[enrichments$Pval <0.05, "ci_right" ], na.rm = T))) 

