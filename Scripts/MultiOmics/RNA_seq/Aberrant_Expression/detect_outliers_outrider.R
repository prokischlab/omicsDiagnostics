#'---
#' title: OUTRIDER pipeline with OHT implementation
#' author: Dmitrii Smirnov
#' wb:
#'  input:
#'   - sa: '`sm config["ANNOTATION"]`'
#'   - gencode_annotation: '`sm config["DATASETS"] + "/gene_annotation_v29.tsv"`'
#'   - counts: '`sm config["RAW_RNA"]`'
#'   - txdb: '`sm config["DATASETS"] + "/txdb.db"`'
#'  output:
#'   - ods_unfiltered: '`sm config["PROC_DATA"] + "/outrider/ods_unfiltered.Rds"`'
#'   - ods: '`sm config["PROC_DATA"] + "/outrider/ods.Rds"`'
#'   - results: '`sm config["PROC_DATA"] + "/outrider/OUTRIDER_results.rds"`'
#'  threads: 8
#' output: 
#'   html_document:
#'    code_folding: hide
#'    code_download: TRUE
#'---

suppressPackageStartupMessages({
  library(ggplot2)
  library(data.table)
  library(dplyr)
  library(magrittr)
  library(readr)
  library(tibble)
  library(SummarizedExperiment)
  library(OUTRIDER)
})

###############################
# Part 1: Preprocessing & Filtering
###############################

fpkmCutoff <- 1

# Load counts data and set gene IDs as rownames
counts <- fread(snakemake@input$counts) %>% as.data.frame()
rownames(counts) <- counts$gene_id
counts$gene_id <- NULL

# Prepare sample annotation for OUTRIDER (subset to selected samples)
sa <- fread(snakemake@input$sa)
sa <- sa[USE_FOR_PROTEOMICS_PAPER == TRUE]
sa <- sa[ SAMPLE_ID %in% colnames(counts), .(SAMPLE_ID, gender)]
sa[, sampleID := SAMPLE_ID]
sa$SAMPLE_ID <- NULL
sa <- as.data.frame(sa)
rownames(sa) <- sa$sampleID

# Add gene annotation (rowData)
gene_annot <- fread(snakemake@input$gencode_annotation)
gene_annot_data <- data.table(gene_id_unique = rownames(counts))
gene_annot_data <- left_join(gene_annot_data,
                             gene_annot[, .(gene_id_unique, gene_name_unique, gene_type, gene_status)],
                             by = "gene_id_unique")

# Create Outrider dataset object
ods <- OutriderDataSet(countData = counts, colData = sa, rowData = gene_annot_data)

# Filter genes using txdb (using UCSC style)
gencode_txdb <- loadDb(snakemake@input$txdb)
seqlevelsStyle(gencode_txdb) <- "UCSC"
gencode_txdb <- keepStandardChromosomes(gencode_txdb)

# Apply filtering (calculate FPKM and mark low expressed genes; filter=FALSE returns full object)
ods <- filterExpression(ods, gtfFile = gencode_txdb, filter = FALSE, fpkmCutoff = fpkmCutoff)
plotFPKM(ods) + theme_bw(base_size = 14)

# Mark genes with at least one nonzero count
rowData(ods)$counted1sample <- rowSums(assay(ods)) > 0

# Save unfiltered dataset (before further filtering)
saveRDS(ods, snakemake@output$ods_unfiltered)

###############################
# Part 2: Run OUTRIDER pipeline using OHT (autoencoder) 
###############################

# Subset to genes passing the filtering step
ods <- ods[mcols(ods)$passedFilter, ]

# Estimate size factors as required by OUTRIDER
ods <- estimateSizeFactors(ods)

# With the new OHT implementation, simply calling OUTRIDER without manually specifying q 
# will trigger automatic estimation of the encoding dimension via OHT.
message("Starting OUTRIDER fitting using OHT...")
# Use SerialParam() or MulticoreParam() as needed; here we use SerialParam() for determinism.
ods <- OUTRIDER(ods, implementation = 'autoencoder', BPPARAM = MulticoreParam(workers=as.integer(snakemake@threads)))
message("OUTRIDER fitting finished.")

# Save the final OUTRIDER object with results
saveRDS(ods, snakemake@output$ods)

# Plot normalized count correlation heatmap after model fitting
plotCountCorHeatmap(ods, normalize = TRUE)
rownames(ods) <- rowData(ods)$gene_name_unique

# Extract results and format output table
res <- results(ods, all = TRUE)
res[, FC := round(2^l2fc, 2)]
res[, geneID := toupper(geneID)]

setnames(res,
         c("geneID", "sampleID", "pValue", "padjust", "normcounts",
           "zScore", "l2fc", "aberrant", "FC"),
         c("geneID", "SAMPLE_ID", "RNA_PVALUE", "RNA_PADJ", "normcounts",
           "RNA_ZSCORE", "RNA_LOG2FC", "RNA_outlier", "RNA_FC"),
         skip_absent = TRUE)

# Correct gene naming if needed
res$geneID[grep('ZNF503_2', res$geneID)] <- 'ZNF503'

# Order and deduplicate the results table
res <- res[, c("SAMPLE_ID", "geneID", "normcounts", "RNA_FC", "RNA_LOG2FC", 
               "RNA_ZSCORE", "RNA_PVALUE", "RNA_PADJ", "RNA_outlier")]
res <- res[!duplicated(res), ]
res <- res[order(res$RNA_PADJ), ]
res <- res[!duplicated(res[, c("SAMPLE_ID", "geneID")]), ]
res <- res[!is.na(res$geneID) & !is.na(res$SAMPLE_ID), ]

# Save the final results table
saveRDS(res, snakemake@output$results)
