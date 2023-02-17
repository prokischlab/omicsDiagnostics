#'---
#' title: Proportion of outliers with rare variants
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





#####################



# As soon as it is not allowed to share genetic data,
# only the result of the analysis could be shared.
# to reproduce the analysis please prepare the variant table (vt)
# with 3 columns: sample_Id, gene_Id and variant type: stop, coding, splice ...



# Read variant annotation
# vt <- fread("../rare_variants.tsv")
# vt <- vt[ geneID %in% unique(rp$geneID), c("SAMPLE_ID", "geneID", "var_type")]
# 
# rp <- merge(vt, rp, by = c("SAMPLE_ID", "geneID"), all.y = T )
# rp[is.na(var_type), var_type := "no rare variant"]
# 
# ####################
# # READ enrichment
# enr <- fread("../enrichment_rare_vt.tsv")
# enr <- enr[, c("var_type", "up_down_outlier", "importance")]
# 
# 
# pr <- merge(rp, enr, by = c("var_type", "up_down_outlier"), all.x=T)
# pr <- pr[order(importance), ]
# pr <- pr[!duplicated(pr[, c("up_down_outlier",  "SAMPLE_ID", "geneID")]), ]
# pr$importance <-NULL
# pr <- pr[!duplicated(pr), ]
# 
# 
# # Calculate proportions of outliers
# pr <- pr[, .N, by = .(up_down_outlier, var_type)]
# pr[, total := sum(N), by = .(up_down_outlier)]
# pr[, prop := N/sum(N), by = up_down_outlier]
# pr[, type := "var_type"]
# write_tsv(pr, "../proportions_rare.tsv")


####################
pr <- fread(snakemake@input$enrichments_proportions)
pr <- pr[ type == "var_type"]

#' ### Underexpression outliers
pr[, var_type := factor(var_type, levels = c("no rare variant", "non_coding", "synonymous", "coding", "frameshift", "splice", "stop"))]
pr[, outlier_class_label := paste0(up_down_outlier ,'\n', "(n = ", total, ")") ]

#+ fig.width=13, fig.height=3
ggplot(pr[up_down_outlier %in% c("RNA_underexpression", "Protein_underexpression", "RNA_Protein_underexpression", "non_outlier")], aes(outlier_class_label, prop)) +
  geom_bar(stat= 'identity', aes(fill = var_type)) +
  scale_fill_manual(values = variant_colors ) +
  coord_flip(ylim = c(0,1)) +
  scale_y_continuous(breaks=seq(0,1, 0.1), labels=scales::percent) +
  labs( y = "Proportion of outliers with rare variants")+
  theme_classic()+
  theme(legend.position="top",  axis.title.y = element_blank() ,legend.direction = "horizontal",
        legend.title = element_blank(),
        axis.text.y = element_text(face="bold", size=12)) +
  guides(fill = guide_legend(nrow = 1))


