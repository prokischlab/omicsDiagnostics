#'---
#' title: Integrate OMICs       
#' author: Dmitrii Smirnov
#' wb:
#'  log:
#'   - snakemake: '`sm config["PROC_DATA"] + "/integration/integration_log.RDS"`'
#'  input:
#'  - config: 'src/config.R'
#'  - sample_annotation: '`sm config["ANNOTATION"]`'
#'  - outrider_results: '`sm config["PROC_DATA"] + "/outrider/OUTRIDER_results.rds"`'
#'  - protrider_results: '`sm config["PROC_DATA"] + "/protrider/PROTRIDER_results.rds"`'
#'  - phenotype_data: '`sm config["PROC_DATA"] + "/HPO/Patients_phenotype_data.tsv"`'
#'  - var_hpo: '`sm config["RAW_DATA"] + "/patient_variant_hpo_data.tsv"`'
#'  output:
#'  - patient_omics_full: '`sm config["PROC_DATA"] + "/integration/patient_omics_full.RDS"`'
#'  - patient_omics: '`sm config["PROC_DATA"] + "/integration/patient_omics.RDS"`'
#' output: 
#'   html_document:
#'    code_folding: hide
#'    code_download: TRUE
#'---

saveRDS(snakemake, snakemake@log$snakemake)

# snakemake <- readRDS("/Users/Mitya/Desktop/working/omicsDagnostics_data/processed_data/integration/integration_log.RDS")

# Load config and functions
source(snakemake@input$config)
source("src/functions/Integration/integrate_annotate_omics.R")


# Load sample annotation
sa <- fread(snakemake@input$sample_annotation)
sa <- sa[USE_FOR_PROTEOMICS_PAPER == T]

# Load OUTRIDER results
rna <- readRDS(snakemake@input$outrider_results) %>% as.data.table()

head(rna)

# Load PROTRIDER results
prot <- readRDS(snakemake@input$protrider_results) %>% as.data.table()
head(prot)

# Load Phenotype data - only for diagnosed cases
phenotype <- fread(snakemake@input$phenotype_data)
head(phenotype)

# Load variant (WES) and phenotype information for all cases
var_hpo <- fread(snakemake@input$var_hpo)
head(var_hpo)

# Combine OMICs
patient_omics <- combine_omics(annotation = sa,
                                  RNA_ae = rna, 
                                  PROTEIN_ae = prot, 
                                  variants_hpo = var_hpo,
                                  Padj_threshold =  PADJ_THRESHOLD,
                                  Zscore_threshold = ZSCORE_THRESHOLD, 
                                  l2FC_threshold = LOG2FC_THRESHOLD )


head(patient_omics)
patient_omics[ , MDG := ""]

# Save full results
saveRDS(patient_omics, snakemake@output$patient_omics_full)

# Select and rename columns
selected_cols <- c("SAMPLE_ID", "geneID", "causal_gene",
                  "rare", "potential_biallelic", "gene_class",
                  "HPO_match", "Semantic_sim", "MDG", 
                  "outlier_class", "validated", "gene_detected",
                  "RNA_ZSCORE", "PROTEIN_ZSCORE",
                  "RNA_FC", "PROTEIN_FC", 
                  "RNA_PVALUE", "PROTEIN_PVALUE",
                  "rank_rna", "normcounts",
                  "rank_protein", "PROTEIN_INT",
                  "WES_avaliable", "HPO_avaliable", 
                  "RNA_seq_avaliable", "Proteomics_avaliable")

# Ensure all columns exist
missing_cols <- setdiff(selected_cols, names(patient_omics))
if (length(missing_cols) > 0) {
  warning("Missing columns: ", paste(missing_cols, collapse = ", "))
  # Add missing columns with NA values
  for (col in missing_cols) {
    patient_omics[[col]] <- NA
  }
}

# Select and reorder columns
patient_omics <- patient_omics[, ..selected_cols]

# Save results
saveRDS(patient_omics, snakemake@output$patient_omics)




#' ## Number of outliers per class
os <- patient_omics[, .N, by = outlier_class]
paste("N individuals = ",  uniqueN(patient_omics$SAMPLE_ID))
#+echo=F
DT::datatable(os, 
              caption = "Number of outliers per class",  escape = F)


range <- c(patient_omics$RNA_ZSCORE, patient_omics$PROTEIN_ZSCORE)
ggplot(patient_omics, aes(RNA_ZSCORE, PROTEIN_ZSCORE)) +
  geom_point( aes( color= outlier_class)) +
  #stat_cor( aes(label = paste(..r.label.. , ..rr.label.., ..p.label.., sep = "~`,`~")),  na.rm = T, method = "spearman")+
  xlab("RNA zScore") + 
  ylab("Protein zScore") +
  geom_vline(xintercept = 0, color = "grey50") +
  geom_hline(yintercept = 0, color = "grey50") +
  geom_abline(color = "grey60")+#
  scale_color_manual(values = outlier_colors)+ 
  coord_fixed( xlim = c( -max(abs(range) , na.rm = T) , max(abs(range) , na.rm = T) ),  ylim = c( -max(abs(range) , na.rm = T),  max(abs(range) , na.rm = T)  )) + 
  ggtitle("OUTRIDER VS PROTRIDER")+
  theme_bw()



#' ## Number of outliers per sample
patient_omics <- patient_omics[ outlier_class != "non_outlier"]

# count outliers by sample
os <- patient_omics[ ,  .N, by = SAMPLE_ID]
paste("Samples with at least 1 outlier = ",  uniqueN(os$SAMPLE_ID))

# Create reference table with all outlier classes
osx <- data.table(
  SAMPLE_ID = rep(os$SAMPLE_ID, each = 3),
  outlier_class = rep(c("RNA-only", "protein-only", "RNA-and-protein"), times = length(os$SAMPLE_ID))
)

# Merge with actual counts
os <- merge(os, osx, by = "SAMPLE_ID", all = TRUE)
setnames(os, c("SAMPLE_ID", "Outliers_per_sample", "outlier_class"))

# Count by outlier class and sample 
osx <- patient_omics[,  .N, by = .(outlier_class, SAMPLE_ID)]
os <- merge(os, osx, by = c("SAMPLE_ID", "outlier_class"), all.x = TRUE)
os[is.na(N), N := 0]
os <- os[order(Outliers_per_sample),]
os[, Rank := 1:.N, by = outlier_class]



ggplot(os, aes(Rank , N, color = outlier_class))+
  geom_line( size = 1)+
  theme_classic()+
  scale_color_manual(values = outlier_colors)+ 

  scale_y_continuous(trans='log10')+

  xlab("Sample rank") + 
  ylab("# aberrantly expressed genes")+
  theme(plot.title = element_text(hjust = 0.5),
        axis.title.x= element_text( size=12, margin = NULL,face="bold"),
        axis.title.y= element_text( size=12, margin = NULL,face="bold"),
        axis.text.x = element_text(face="bold",  size=12),
        axis.text.y = element_text(face="bold",  size=12),
        legend.title = element_blank(),
        legend.position = c(0.2, 0.8))




#' ## Results table
#+echo=F
DT::datatable(patient_omics, 
              caption = "Outrider VS Protrider results", style = 'bootstrap', filter = 'top', escape = F,
              extensions = c( 'Buttons', 'ColReorder' ), 
              options = list( colReorder = TRUE, dom = 'Bfrtip',
                              buttons = c('copy', 'csv', 'excel', 'pdf', 'print')))



