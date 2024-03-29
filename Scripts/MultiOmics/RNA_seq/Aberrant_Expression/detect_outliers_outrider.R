#'---
#' title: OUTRIDER pipeline
#' author: Michaela Mueller, Vicente Yepez, Dmitrii Smirnov
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
#'  threads: 40
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



#' # Part1: Filter counts
#########################################################
fpkmCutoff <-  1


# Get counts data
counts <- fread(snakemake@input$counts) %>%  as.data.frame()
rownames(counts) <- counts$gene_id
counts$gene_id <- NULL



# Add sample annotation for OUTRIER object colData
sa <- fread(snakemake@input$sa)
sa <- sa[USE_FOR_PROTEOMICS_PAPER == T]
sa <- sa[ SAMPLE_ID %in% colnames(counts) , .(SAMPLE_ID, gender)]
sa[, sampleID := SAMPLE_ID]
sa$SAMPLE_ID <- NULL
sa <- as.data.frame(sa)
rownames(sa) <- sa$sampleID


# Add gene annotation rowData
gene_annot <- fread(snakemake@input$gencode_annotation)
gene_annot_data <- data.table(gene_id_unique = rownames(counts))
gene_annot_data <- left_join(gene_annot_data, gene_annot[,.(gene_id_unique, gene_name_unique, gene_type, gene_status)], by = "gene_id_unique")

# Create Outrider dataset obj
ods <- OutriderDataSet(countData=counts, colData = sa, rowData = gene_annot_data)


# filter not expressed genes
gencode_txdb <- loadDb(snakemake@input$txdb)
seqlevelsStyle(gencode_txdb) <- "UCSC"
gencode_txdb <- keepStandardChromosomes(gencode_txdb)

ods <- OUTRIDER::filterExpression(ods, gtfFile=gencode_txdb, filter=FALSE, fpkmCutoff=fpkmCutoff)
g <- plotFPKM(ods) + theme_bw(base_size = 14)

#' ## Count filtering results, fpkm cutoff = 1
g

rowData(ods)$counted1sample = rowSums(assay(ods)) > 0


# Save the ods object before filtering, so as to preserve the original number of genes
saveRDS(ods, snakemake@output$ods_unfiltered)



#' # Part 2: RUN OUTRIDER pipeline
#########################################################
# threads <-  30

ods <- ods[mcols(ods)$passedFilter,] 
ods <- OUTRIDER::estimateSizeFactors(ods)

a = 5 
b = min(ncol(ods), nrow(ods)) / 3   # N/3
Nsteps = min(20, ncol(ods)/3, nrow(ods)/3)   # Do either 20 steps or N
pars_q <- round(exp(seq(log(a),log(b),length.out = Nsteps))) %>% unique 



#' ## Before OUTRIDER
plotCountCorHeatmap(ods, normalize=F)

message("outrider fitting start")
register(MulticoreParam(workers=as.integer(snakemake@threads)))
ods <- OUTRIDER::findEncodingDim(ods, lnorm = T, params = pars_q, implementation = 'autoencoder', BPPARAM = SerialParam() ) # BPPARAM = MulticoreParam(workers=as.integer(snakemake@threads)) , BPPARAM = SerialParam()
plotEncDimSearch(ods)
message("encoding dimention found")
ods <- OUTRIDER::OUTRIDER(ods, implementation = 'autoencoder', BPPARAM = SerialParam()) #  ,  BPPARAM = MulticoreParam(workers=as.integer(snakemake@threads))
message("outrider fitting finished")


# Save  OUTRIDER OBJECT WITH RESULTS
saveRDS(ods, snakemake@output$ods)



#' ## After OUTRIDER
plotCountCorHeatmap(ods, normalize=TRUE)
row.names(ods) <- rowData(ods)$gene_name_unique


# Save  OUTRIDER OBJECT WITH RESULTS
saveRDS(ods, snakemake@output$ods)


res <- OUTRIDER::results(ods, all = TRUE)
res[, FC := round(2^l2fc, 2)]
res[, geneID := toupper(geneID)]


setnames(res,
         c("geneID", "sampleID", "pValue", "padjust", "normcounts",
           "zScore",   "l2fc",    "aberrant", "FC"), 
         c("geneID" , "SAMPLE_ID", "RNA_PVALUE" , "RNA_PADJ", "normcounts",
           "RNA_ZSCORE","RNA_LOG2FC", "RNA_outlier",  "RNA_FC" ) , skip_absent=TRUE) 
res$geneID[grep('ZNF503_2', res$geneID)] <- 'ZNF503'  # This gene appears as ZNF503_2



res <- res[,c("SAMPLE_ID", "geneID", "normcounts", "RNA_FC","RNA_LOG2FC", "RNA_ZSCORE", "RNA_PVALUE" , "RNA_PADJ","RNA_outlier")] 
res <- res[!duplicated(res), ]
res <- res[order(res$RNA_PADJ), ]
res <-res[!duplicated(res[, c("SAMPLE_ID", "geneID")]),  ]

res <- res[!is.na(res$geneID), ]
res <- res[!is.na(res$SAMPLE_ID), ]


saveRDS(res,  snakemake@output$results)







