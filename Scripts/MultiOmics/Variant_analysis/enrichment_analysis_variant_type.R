#'---
#' title: Enrichment of rare variant types in outliers  
#' author: smirnovd
#' wb:
#'  input:
#'  - patient_omics: '`sm config["PROC_DATA"] + "/integration/patient_omics_full.RDS"`'
#'  - enrichments_proportions: '`sm config["RAW_DATA"] + "/enrichment_proportions_variants.tsv"`'
#'  output:
#' output: 
#'   html_document:
#'    code_folding: hide
#'    code_download: TRUE
#'---



# Load plotting config and functions
source("src/config.R")
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

# Subset necessary columns
rp <- rp[ , c("SAMPLE_ID", "geneID", "outlier_class", "up_down_outlier") ]
rp <- rp[!duplicated(rp), ]

#' # Outlier class
os1 <- rp[, .N, by = outlier_class]
os2 <- rp[, .N, by = up_down_outlier]
colnames(os2)[1] <-"outlier_class" 
os <- rbind(os1, os2)
rm(os1, os2)
#' ### Number of outliers per class
DT::datatable(os, caption = "Outlier counts",   escape = F)



####################


# As soon as it is not allowed to share genetic data,
# only the result of the analysis could be shared.
# to reproduce the analysis please prepare the variant table (vt)
# with 3 columns: sample_Id, gene_Id and variant type: stop, coding, splice ...

# Read variant annotation 
# vt <- fread("../rare_variants.tsv")
# vt <- vt[ geneID %in% unique(rp$geneID), c("SAMPLE_ID", "geneID", "var_type")]
# rp <- rp[SAMPLE_ID %in% unique(vt$SAMPLE_ID) ]
# rp <- merge(vt, rp, by = c("SAMPLE_ID", "geneID"), all.y = T )
# rp[is.na(var_type), var_type := "no rare variant"]
# 
# 
# enr <- reshape2::dcast(data = rp, geneID +  SAMPLE_ID +  var_type ~ up_down_outlier , fun.aggregate = cat.agg)
# 
# enrichments <- data.frame()
# for (outl in unique(rp$up_down_outlier)){
#   tem_df <- enr
#   setnames(tem_df, outl, "outliers")
#   tem_df <- tem_df[, c("geneID", "SAMPLE_ID", "outliers", "var_type" )]
#   feat <- unique(tem_df$var_type)
#   tem_df <- dcast(data = tem_df, geneID + SAMPLE_ID  + outliers ~ var_type, fun.aggregate = cat.agg) 
#   enrichment <- enrich(tem_df, outlier= "outliers" ,features= feat)
#   enrichment$ci_left <- enrichment$Estim - 1.96 * enrichment$Std
#   enrichment$ci_right <- enrichment$Estim + 1.96 * enrichment$Std
#   enrichment$outlier_class <- rep(outl, nrow(enrichment))
#   enrichments <- rbind(enrichments, enrichment)
#   enr$outliers <- NULL
# }
# 
# enrichments <- as.data.table(enrichments)
# enrichments[, significant := 1 ]
# enrichments[enrichments$Pval >= 0.05 , significant := 0.9 ]
# 
# 
# enrichments <- as.data.table(enrichments)
# enrichments[outlier_class == "non_outlier", priori:=Estim]
# enrichments[ outlier_class != "non_outlier" & Pval >= 0.05, priori:=Estim]
# enrichments[ outlier_class != "non_outlier" & Pval < 0.05, priori:=Estim +1]
# enrichments <- enrichments[order(priori, decreasing = T), ]
# enrichments$importance <- seq(1:nrow(enrichments))
# enrichments <- enrichments[order(importance), ]
# 
# setnames(enrichments, c("Cat", "outlier_class"), c("var_type", "up_down_outlier"))
# write_tsv(enrichments, "../enrichment_rare_vt.tsv")
# 


enrichments <- fread(snakemake@input$enrichments_proportions)
# enrichments <- fread("/s/project/mitoMultiOmics/multiOMICs_integration/raw_data/enrichment_proportions_variants.tsv")
enrichments <- enrichments[ type == "var_type"]
enrichments[, var_type := factor(var_type, levels = c("no rare variant", 
                                             "non_coding", "synonymous", "coding", "frameshift", "splice", "stop"))]

nodelist <- levels(enrichments$var_type)



#' ## Outlier class
#+ fig.width=6, fig.height=9
ggplot(data = enrichments[up_down_outlier %in% c("RNA_underexpression", "Protein_underexpression", "RNA_Protein_underexpression", "non_outlier")], 
       aes(x = var_type, y = Estim)) +
  geom_pointrange(aes(x = var_type, ymin = ci_left, ymax = ci_right, colour = up_down_outlier, alpha = significant), position = position_dodge(width = 1)) +
  geom_hline(yintercept = 0) +
  geom_vline(color = "gray80", linetype = "dashed", xintercept = which(nodelist=="non_coding")-0.5) +
  geom_vline(color = "gray80", linetype = "dashed", xintercept = which(nodelist=="synonymous")-0.5) +
  geom_vline(color = "gray80", linetype = "dashed", xintercept = which(nodelist=="coding")-0.5) +
  geom_vline(color = "gray80", linetype = "dashed", xintercept = which(nodelist=="frameshift")-0.5) +
  geom_vline(color = "gray80", linetype = "dashed", xintercept = which(nodelist=="splice")-0.5) +
  geom_vline(color = "gray80", linetype = "dashed", xintercept = which(nodelist=="stop")-0.5) +
  theme_bw() + 
  xlab('') + 
  ylab('Log odds ratio')  + 
  scale_alpha(guide = 'none')+
  scale_color_manual(values = outlier_colors)+ 
  coord_flip(ylim = c(min(enrichments[enrichments$Pval <0.05, 'ci_left' ] , na.rm = T),  max(enrichments[enrichments$Pval <0.05, "ci_right" ], na.rm = T))) 



