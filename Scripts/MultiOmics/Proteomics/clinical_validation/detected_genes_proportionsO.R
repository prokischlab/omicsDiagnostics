#'---
#' title: Gene detection by OMICs proportions
#' author: smirnovd
#' wb:
#'  input:
#'  - patient_omics: '`sm config["PROC_DATA"] + "/integration/patient_omics.RDS"`'
#'  - gencode_annotation: '`sm config["DATASETS"] + "/gene_annotation_v29.tsv"`'
#'  - detected_transcripts: '`sm config["PROC_DATA"] + "/integration/detected_transcripts.csv"`'
#'  - detected_proteins: '`sm config["PROC_DATA"] + "/integration/detected_proteins.csv"`'
#'  - detected_proteins_gtex: '`sm config["PROC_DATA"] + "/integration/detected_proteins_gtex.csv"`'
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


allgenes_annotation <- merge(genecode, dg, by = "geneID", all.x = T   )
rm(genecode, dg)


# Load list of genes, detected by RNS-seq
# detected_transcripts <- fread('/s/project/mitoMultiOmics/multiOMICs_integration/processed_data/integration/detected_transcripts.csv')
detected_transcripts <- fread(snakemake@input$detected_transcripts)
detected_transcripts[ , RNA := det_rate]
detected_transcripts <- detected_transcripts[ , .(geneID, RNA)]
detected_transcripts <- detected_transcripts[!duplicated(detected_transcripts)]
detected_transcripts[ , TISSUE := "Skin fibroblasts"]
detected_transcripts[ , Source := "GENOMIT"]



# Load list of genes, detected by proteomics
# detected_proteins <- fread('/s/project/mitoMultiOmics/multiOMICs_integration/processed_data/integration/detected_proteins.csv')
detected_proteins <- fread(snakemake@input$detected_proteins)
detected_proteins[ , Protein := det_rate]
detected_proteins <- detected_proteins[ , .(geneID, Protein)]
detected_proteins <- detected_proteins[!duplicated(detected_proteins)]
detected_proteins[ , TISSUE := "Skin fibroblasts"]
detected_proteins[ , Source := "GENOMIT"]

# Load list of genes, detected by GTEx proteomics
# detected_proteins_gtex <- fread('/s/project/mitoMultiOmics/multiOMICs_integration/processed_data/integration/detected_proteins_gtex.csv')
detected_proteins_gtex <- fread(snakemake@input$detected_proteins_gtex)
detected_proteins_gtex[ , Protein := det_rate]
detected_proteins_gtex <- detected_proteins_gtex[ , .(geneID, Protein, TISSUE)]
detected_proteins_gtex <- detected_proteins_gtex[!duplicated(detected_proteins_gtex)]
detected_proteins_gtex <- detected_proteins_gtex[!is.na(Protein)]
detected_proteins_gtex[ , Source := "GTEx"]

detected_proteins_tissue <- rbind(detected_proteins,detected_proteins_gtex)

detected_omics <- merge(detected_transcripts, detected_proteins_tissue, by = c("geneID", "TISSUE", "Source"), all = T) 
#

allgenes <- as.data.table( merge(allgenes_annotation, detected_omics, by = "geneID", all.x = T) )
allgenes <- allgenes[!is.na( TISSUE )]
allgenes[is.na(Protein), Protein := "not detected" ]

allgenes[is.na(RNA) & TISSUE == "Skin fibroblasts", RNA := "not detected" ]

allgenes <- allgenes[order( Source,  disease)]

write_csv(allgenes, '/s/project/mitoMultiOmics/multiOMICs_integration/processed_data/integration/detected_omics.csv')





#+echo=F
DT::datatable(allgenes, caption = "Gene coverage by omics", 
              style = 'bootstrap', filter = 'top', escape = F,
              extensions = c( 'Buttons', 'ColReorder' ), 
              options = list( colReorder = TRUE, dom = 'Bfrtip',
                              buttons = c('copy', 'csv', 'excel', 'pdf', 'print')))


detected_proteins_gtex[ , gtex_protein := paste0( TISSUE, ": ",  Protein)]
detected_proteins_gtex <- detected_proteins_gtex[ , .(geneID, gtex_protein)]
detected_proteins_gtex <- detected_proteins_gtex[!duplicated(detected_proteins_gtex)]
detected_proteins_gtex[ , aggregated_tissues := toString(gtex_protein), by = geneID]


dp_gtex <- detected_proteins_gtex[ , c("geneID" , "aggregated_tissues")]
dp_gtex <- dp_gtex[!duplicated(dp_gtex)]
setnames(dp_gtex, c("geneID" , "gtex_protein") )



detected_omics2 <- merge(detected_transcripts, detected_proteins, by = c("geneID", "TISSUE", "Source"), all = T) 
detected_omics2 <- merge(detected_omics2, dp_gtex, by = c("geneID"), all = T) 
detected_omics2$TISSUE <- NULL
detected_omics2$Source <- NULL
detected_omics2[is.na(Protein), Protein := "not detected" ]
detected_omics2[is.na(RNA), RNA := "not detected" ]
detected_omics2[is.na(gtex_protein), gtex_protein := "not detected" ]

allgenes2 <- as.data.table( merge(allgenes_annotation, detected_omics2, by = "geneID") )
allgenes2 <- allgenes2[order(  disease, Protein)]

write_csv(allgenes2, 'S_Table4_detected_omics.csv')





