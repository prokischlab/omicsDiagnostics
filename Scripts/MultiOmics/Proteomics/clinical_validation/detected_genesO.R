#'---
#' title: Gene detection by OMICs
#' author: smirnovd
#' wb:
#'  input:
#'  - patient_omics: '`sm config["PROC_DATA"] + "/integration/patient_omics.RDS"`'
#'  - gencode_annotation: '`sm config["DATASETS"] + "/gene_annotation_v29.tsv"`'
#'  - detected_transcripts: '`sm config["PROC_DATA"] + "/integration/detected_transcripts.tsv"`'
#'  - detected_proteins: '`sm config["PROC_DATA"] + "/integration/detected_proteins.tsv"`'
#'  - detected_proteins_gtex: '`sm config["PROC_DATA"] + "/integration/detected_proteins_gtex.tsv"`'
#'  - disease_genes: '`sm config["DATASETS"] + "/disease_genes.tsv"`'
#'  output:
#' output: 
#'   html_document:
#'    code_folding: hide
#'    code_download: TRUE
#'---

# Load config
source("src/config.R")


# Read integrated omics file 
# rp <- readRDS("/s/project/mitoMultiOmics/multiOMICs_integration/processed_data/integration/patient_omics.RDS") %>% as.data.table()
rp <- readRDS(snakemake@input$patient_omics) %>% as.data.table()

rp <- rp[ gene_detected != "no RNA"]
rp <- rp[ !is.na(normcounts )]

ggplot(rp, aes(log10(normcounts), fill = gene_detected))+
  geom_density(alpha = 0.4)+
  theme_bw()+
  ggtitle("RNA counts for non-detected proteins")


rm(rp)

# Get all protein coding genes
# genecode_v29 <- fread('/s/project/mitoMultiOmics/multiOMICs_integration/datasets/gene_annotation_v29.tsv')
genecode_v29 <- fread(snakemake@input$gencode_annotation)
genecode_v29[, geneID := toupper(gene_name)]
genecode_v29 <- genecode_v29[ , .(geneID , gene_type)]
genecode_v29 <- genecode_v29[!duplicated(genecode_v29)]
genecode_v29 <- genecode_v29[order(gene_type)]

genecode_v29[ , gencode_v29 := toString(gene_type), by = geneID]
genecode <- genecode_v29[, c("geneID" , "gencode_v29")]
genecode <- genecode[!duplicated(genecode), ]
rm(genecode_v29)

# Load disease genes table
# dis_genes <- fread('/s/project/mitoMultiOmics/multiOMICs_integration/datasets/disease_genes.tsv')
dis_genes <- fread(snakemake@input$disease_genes)
dis_genes <- dis_genes[ , .(geneID , DISEASE)]
dis_genes <- dis_genes[!duplicated(dis_genes)]
dis_genes <- dis_genes[order(DISEASE)]

dis_genes[ , disease := toString(DISEASE), by = geneID]
dg <- dis_genes[ , c("geneID" , "disease")]
dg <- dg[!duplicated(dg), ]
rm(dis_genes)


allgenes <- merge(genecode, dg, by = "geneID", all.x = T   )
rm(genecode, dg)


# Load list of genes, detected by RNS-seq
# detected_transcripts <- fread('/s/project/mitoMultiOmics/multiOMICs_integration/processed_data/integration/detected_transcripts2.tsv')
detected_transcripts <- fread(snakemake@input$detected_transcripts)
detected_transcripts[once == T , fib_RNA := "once"]
detected_transcripts[half == T , fib_RNA := "half of the samples"]
detected_transcripts[all == T , fib_RNA := "all of the samples"]
detected_transcripts <- detected_transcripts[ , .(geneID, fib_RNA)]
detected_transcripts <- detected_transcripts[!duplicated(detected_transcripts)]


allgenes <- as.data.table( merge(allgenes, detected_transcripts, by = "geneID", all.x = T) )
allgenes[is.na(fib_RNA), fib_RNA := "not detected" ]

# Load list of genes, detected by proteomics
# detected_proteins <- fread('/s/project/mitoMultiOmics/multiOMICs_integration/processed_data/integration/detected_proteins.tsv')
detected_proteins <- fread(snakemake@input$detected_proteins)
detected_proteins[once == T , fib_protein := "once"]
detected_proteins[half == T , fib_protein := "half of the samples"]
detected_proteins[all == T , fib_protein := "all of the samples"]
detected_proteins <- detected_proteins[ , .(geneID, fib_protein)]
detected_proteins <- detected_proteins[!duplicated(detected_proteins)]
allgenes <- as.data.table( merge(allgenes, detected_proteins, by = "geneID", all.x = T) )
allgenes[is.na(fib_protein), fib_protein := "not detected" ]


# Load list of genes, detected by GTEx proteomics
# detected_proteins_gtex <- fread('/s/project/mitoMultiOmics/multiOMICs_integration/processed_data/integration/detected_proteins_gtex.tsv')
detected_proteins_gtex <- fread(snakemake@input$detected_proteins_gtex)
detected_proteins_gtex[once == T , protein := "once"]
detected_proteins_gtex[half == T , protein := "half of the samples"]
detected_proteins_gtex[all == T , protein := "all of the samples"]
detected_proteins_gtex <- detected_proteins_gtex[ , .(geneID, protein, TISSUE)]
detected_proteins_gtex <- detected_proteins_gtex[!is.na(protein)]
detected_proteins_gtex[ , gtex_protein := toString(TISSUE), by = .(geneID, protein)]
detected_proteins_gtex$TISSUE <- NULL
detected_proteins_gtex <- detected_proteins_gtex[!duplicated(detected_proteins_gtex)]

detected_proteins_gtex[ , gtex_protein := paste(protein , "in",  gtex_protein)]

detected_proteins_gtex <- detected_proteins_gtex[ , .(geneID, gtex_protein)]
detected_proteins_gtex <- detected_proteins_gtex[!duplicated(detected_proteins_gtex)]

detected_proteins_gtex[ , aggregated_tissues := toString(gtex_protein), by = geneID]
dp_gtex <- detected_proteins_gtex[ , c("geneID" , "aggregated_tissues")]
dp_gtex <- dp_gtex[!duplicated(dp_gtex)]
setnames(dp_gtex, c("geneID" , "gtex_protein") )
rm(detected_proteins_gtex)


allgenes <- as.data.table( merge(allgenes, dp_gtex, by = "geneID", all.x = T) )
allgenes[is.na(gtex_protein), gtex_protein := "not detected" ]


allgenes <- allgenes[ !(fib_RNA == "not detected" &  fib_protein == "not detected" & gtex_protein == "not detected" )]

allgenes <- allgenes[order(disease)]

#+echo=F
DT::datatable(allgenes, caption = "Gene coverage by omics", 
              style = 'bootstrap', filter = 'top', escape = F,
              extensions = c( 'Buttons', 'ColReorder' ), 
              options = list( colReorder = TRUE, dom = 'Bfrtip',
                              buttons = c('copy', 'csv', 'excel', 'pdf', 'print')))










