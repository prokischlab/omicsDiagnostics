

hpo <- get_ontology("datasets/hp.obo", extract_tags="everything")


# only solved 
pat_hpo <- fread("../omicsDagnostics_data/raw_data/Patient_HPO_phenotypes.tsv")
pat_hpo[HPO_ID == "HP:0040168", HPO_ID := "HP:0007359" ]
pat_hpo[HPO_ID == "HP:0001163", HPO_ID := "HP:0005916" ]
pat_hpo[HPO_ID == "HP:0001145", HPO_ID := "HP:0000532" ]
pat_hpo[HPO_ID == "HP:0001146", HPO_ID := "HP:0000580" ]
pat_hpo[HPO_ID == "HP:0001150", HPO_ID := "HP:0000610" ]
write_tsv(pat_hpo, "../omicsDagnostics_data/raw_data/Patient_HPO_phenotypes.tsv")





### Extract minimal set 
colnames(pat_hpo)
pat_hpo <- pat_hpo[ , c("SAMPLE_ID", "HPO_ID"  )]
pat_hpo <- pat_hpo[ !duplicated( pat_hpo)]
uniqueN(pat_hpo$SAMPLE_ID)

pat_hpoM <- c()
for ( pat in unique(pat_hpo$SAMPLE_ID)){
  
  hpOfInterest <- pat_hpo[SAMPLE_ID == pat]$HPO_ID
  
  # extract minimal set
  res <- as.data.frame(minimal_set( hpo, hpOfInterest ) )
  colnames(res)[1] <- "HPO_ID"
  res$SAMPLE_ID <- rep( pat, nrow(res))
  pat_hpoM <- rbind(pat_hpoM, res)
}

pat_hpoX <- as.data.table( pat_hpoM)
pat_hpoZ <- pat_hpoX[ , .N ,  by = SAMPLE_ID]


pat_hpoY <- pat_hpoX[SAMPLE_ID %in% pat_hpoZ[N < 3]$SAMPLE_ID ]
pat_hpoX <- pat_hpoX[SAMPLE_ID %in% pat_hpoZ[N >= 3]$SAMPLE_ID ]

# unique(pat_hpoY$HPO_ID)


# get ancestors for cases with less than 3 hpos 
pat_hpoM <- c()
# pat <- "133876"
for ( pat in unique(pat_hpoY$SAMPLE_ID)){
  
  hpOfInterest <- pat_hpoY[SAMPLE_ID == pat]$HPO_ID
  
  # extract ancestors
  res <- as.data.frame(get_ancestors( hpo, hpOfInterest ) )
  colnames(res)[1] <- "HPO_ID"
  res$SAMPLE_ID <- rep( pat, nrow(res))
  pat_hpoM <- rbind(pat_hpoM, res)
}

pat_hpoX <- rbind(pat_hpoX, pat_hpoM)
pat_hpoX <- pat_hpoX[!duplicated(pat_hpoX)]
pat_hpoZ <- pat_hpoX[ , .N ,  by = SAMPLE_ID]


# pat_hpoX <- merge(pat_hpoM, pat_hpo_sa, by = "SAMPLE_ID")

pat_hpoX <- pat_hpoX[SAMPLE_ID %in% unique(pat_hpoZ[ N > 2 ]$SAMPLE_ID) ]
pat_hpoX <-pat_hpoX[!duplicated(pat_hpoX )]
write_tsv(pat_hpoX, "../omicsDagnostics_data/raw_data/Patient_HPO_phenotypes.tsv")









 pat_hpo <- fread("../omicsDagnostics_data/raw_data/Patient_HPO_phenotypes.tsv")

 




#load gene - hpo 
hpo_gene <- fread("datasets/phenotype_to_genes.txt")
# hpo_gene[, ensembl_gene_id := mapIds(org.Hs.eg.db,
#                                      keys = as.character(ncbi_gene_id),
#                                      column = "ENSEMBL",
#                                      keytype = "ENTREZID",
#                                      multiVals = "first")]
# fwrite(hpo_gene,  "datasets/phenotype_to_genes.txt")


setnames(hpo_gene, c("hpo_id", "hpo_name", "ncbi_gene_id", "gene_symbol", "disease_id", "ensembl_gene_id"),
         c("HPO_ID", "HPO_name", "ncbi_gene_id", "geneID", "disease_id", "ensembl_gene_id"))





