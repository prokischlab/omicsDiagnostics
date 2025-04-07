#'---
#' title: HPO gene table with sample phenotypes
#' author: Dmitrii Smirnov
#' wb:
#'  input:
#'  - patient_hpo: '`sm config["RAW_HPO"]`'
#'  output:
#'  - hpo_gene: '`sm config["RAW_DATA"] + "/hpo_gene.txt"`'
#'  type: script
#'---


######################################
###### Update HPO-gene table #########
######################################

source("src/config.R")

# Load patient's HPO 
pat_hpo <- fread(snakemake@input$patient_hpo)
# pat_hpo <- fread("../omicsDagnostics_data/raw_data/Patient_HPO_phenotypes.tsv")

pat_hpo <- pat_hpo[, c("SAMPLE_ID", "HPO_ID")]
pat_hpo <-pat_hpo[!duplicated(pat_hpo )]
setnames(pat_hpo,  "SAMPLE_ID", "geneID")
pat_hpo[ , type := "patient"]


# Load HPO gene from hpo obo

hpo_gene <- fread("datasets/phenotype_to_genes.txt")

# hpo_gene <- fread("http://purl.obolibrary.org/obo/hp/hpoa/phenotype_to_genes.txt")

# hpo_gene[, ensembl_gene_id := mapIds(org.Hs.eg.db,
#                                      keys = as.character(ncbi_gene_id),
#                                      column = "ENSEMBL",
#                                      keytype = "ENTREZID",
#                                      multiVals = "first")]
# fwrite(hpo_gene,  "datasets/phenotype_to_genes.txt")

setnames(hpo_gene, c("hpo_id", "hpo_name", "ncbi_gene_id", "gene_symbol", "disease_id", "ensembl_gene_id"),
         c("HPO_ID", "HPO_name", "ncbi_gene_id", "geneID", "disease_id", "ensembl_gene_id"))

hpo_gene <- hpo_gene[, c("HPO_ID","geneID")]
hpo_gene <- hpo_gene[!duplicated(hpo_gene)]
hpo_gene[ , type := "gene"]

hpo_gene <- rbind(pat_hpo, hpo_gene)
hpo_gene <- hpo_gene[!duplicated(hpo_gene)]


# Add missing cases and genes wo annotation
hpo_gene <- rbind(hpo_gene, data.table( geneID = "TXNIP", type = "gene", HPO_ID = c("HP:0004902", "HP:0001943", "HP:0003355", "HP:0006568", "HP:0003658", "HP:0004322", "HP:0001508","HP:0001252", "HP:0011968")  ) )
hpo_gene <- rbind(hpo_gene, data.table( geneID = "YARS", type = "gene", HPO_ID = c("HP:0001510", "HP:0001263", "HP:0007305", "HP:0000407", "HP:0000639", "HP:0002611", "HP:0001738", "HP:0001943", "HP:0001903", "HP:0000093", "HP:0100806", "HP:0006528")  ) )

hpo_geneX <- hpo_gene[geneID %in% c("MICOS13", "ADPRS") ]

hpo_geneX[geneID == "MICOS13", geneID := "C19ORF70" ]
hpo_geneX[geneID == "ADPRS", geneID := "ADPRHL2" ]
hpo_geneX[ , type := "gene"]

hpo_gene <- rbind(hpo_gene, hpo_geneX)
hpo_gene <- hpo_gene[!duplicated(hpo_gene)]
rm(pat_hpo,  hpo_geneX)

fwrite(hpo_gene,  snakemake@output$hpo_gene)
# fwrite(hpo_gene, "../omicsDagnostics_data/raw_data/hpo_gene.txt")

#  
# hpo_geneA <- c()
# for ( gn in unique(hpo_gene$geneID)){
# 
#   hpOfInterest <- hpo_gene[geneID == gn]$HPO_ID
# 
#   # extract ancestors
#   res <- as.data.frame(get_ancestors( hpo, hpOfInterest ) )
#   colnames(res)[1] <- "HPO_ID"
#   res$geneID <- rep( gn, nrow(res))
#   hpo_geneA <- rbind(hpo_geneA, res)
# }
# 
# hpo_geneA <- hpo_geneA[!duplicated(hpo_geneA) , ]
# 
# fwrite(hpo_geneA, "../omicsDagnostics_data/raw_data/hpo_gene_ancestors.txt")






