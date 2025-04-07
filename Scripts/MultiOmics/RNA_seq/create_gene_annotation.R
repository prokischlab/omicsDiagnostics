#'---
#' title: Create gene annotation gencode v29
#' author: Vicente Yepez, Dmitrii Smirnov
#' wb:
#'  input:
#'  output:
#'   - txdb: '`sm config["DATASETS"] + "/txdb.db"`'
#'   - gencode_annotation: '`sm config["DATASETS"] + "/gene_annotation_v29.tsv"`'
#'  type: script
#'---

# Get path to the datasets folder
datasets_dir <- yaml::read_yaml("wbuild.yaml")$DATASETS

suppressPackageStartupMessages({
  library(readr)
  library(BSgenome.Hsapiens.UCSC.hg19)
  library(GenomicFeatures)
  library(GenomicRanges)
  library(data.table)
  library(magrittr)
  library(dplyr)
  library(tidyr)
})

# Function to obtain number of junctions of each gene
obtain_junctions <- function(introns_gr){
  it_dt <- as.data.table(introns_gr)
  max_rep_genes <- max(sapply(it_dt$gene_id, length))
  is <- separate(it_dt, gene_id, into = paste0("g", 1:max_rep_genes), sep = ",")
  is <- as.data.table(is)
  jc <- melt(is, measure.vars = paste0("g", 1:max_rep_genes), value.name = 'gene_id')
  jc <- as.data.table(jc)
  jc <- jc[!is.na(gene_id)]
  jc[, gene_id := gsub('c(', '', gene_id, fixed = T)]
  jc[, gene_id := gsub('\"', '', gene_id, fixed = T)]
  jc[, gene_id := gsub(' ', '', gene_id, fixed = T)]
  jc[, gene_id := gsub(')', '', gene_id, fixed = T)]
  jc[, gene_id := gsub('\n', '', gene_id, fixed = T)]
  
  jc[, variable := NULL]
  
  junctions_dt <- jc[, .(N_junctions = .N), by = gene_id]
  return(junctions_dt)
}





# Download gencode gtf
# download.file(destfile= paste0( datasets_dir, "/gencode.v29lift37.annotation.gtf.gz" ),  #snakemake@input$datasets
#               url= "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_29/GRCh37_mapping/gencode.v29lift37.annotation.gtf.gz")


## Create txdb object for OUTRIDER
gencode_txdb = makeTxDbFromGFF(paste0(datasets_dir, "/gencode.v29lift37.annotation.gtf.gz" ), format='gtf') # snakemake@input$datasets



# Subset to include only canonical chromosomes
gencode_txdb <- keepStandardChromosomes(gencode_txdb)

# Save txdb
saveDb(gencode_txdb, snakemake@output$txdb)



## Make annotation table with gene names

gtf_dt <- rtracklayer::import(paste0(datasets_dir , "/gencode.v29lift37.annotation.gtf.gz") ) %>% as.data.table # snakemake@input$datasets
gtf_dt <- as.data.table(gtf_dt)

gtf_dt <- gtf_dt[type == "gene", .(seqnames, start, end, strand, gene_id, gene_name, gene_type, gene_status)]
gtf_dt <- gtf_dt[seqnames %in% GenomeInfoDb::standardChromosomes(BSgenome.Hsapiens.UCSC.hg19)]
setnames(gtf_dt, "gene_id", "gene_id_unique")
gtf_dt <- separate(gtf_dt, "gene_id_unique", into = "gene_id", sep = "\\.", remove = F)
 

# Get genes that appear at least twice
dup_genes <- gtf_dt[duplicated(gtf_dt$gene_name), ]$gene_name 

gtf_dt <- as.data.table(gtf_dt)
# Get genes that appear more than twice
repeated_genes <- names(table(gtf_dt[gene_name %in% dup_genes, ]$gene_name)[table(gtf_dt[gene_name %in% dup_genes, gene_name]) > 1])


# Rename duplicate gene names
gtf_dt[ , Ng := 1:.N, by = gene_name] # warning message
gtf_dt[, gene_name_unique := gene_name]
gtf_dt[Ng > 1, gene_name_unique := paste(gene_name, Ng, sep = '_')]
gtf_dt[, Ng := NULL]


# Add number of junctions
introns_gencode <- intronicParts(gencode_txdb, linked.to.single.gene.only = FALSE)
introns_gencode <- as.data.table(introns_gencode)
gencode_junctions <- obtain_junctions(introns_gencode)

gencode_junctions <- as.data.table(gencode_junctions)

gtf_dt <- left_join(gtf_dt, gencode_junctions, by = c('gene_id_unique' = 'gene_id')) %>% as.data.table()
gtf_dt[is.na(N_junctions), N_junctions := 0]


# Add mtDNA column
gtf_dt[, mtDNA := seqnames == "chrM"]


# Save annotation table
write_tsv(gtf_dt, snakemake@output$gencode_annotation)

