# Recon : Laura data, iem codes


mitog <- read_excel("/Users/Mitya/Downloads/Table_2_The Dimensions of Primary Mitochondrial Disorders.XLSX" ) %>% as.data.table
colnames(mitog)
mitog <- mitog[ , c("Ensembl ID", "Gene name", "Category", "Function",
                    "OMIM_phenotype", "OMIM_MIM number")]
setnames(mitog, c("Ensembl ID", "Gene name", "OMIM_MIM number" ),
         c("Ensembl_ID", "Gene_name", "OMIM_MIM_number" ))
mitog <- mitog[!duplicated(mitog)]
mitog <- mitog[!duplicated(Gene_name)]


ma <- fread("../working/annotations_master/mitoNET_samp_annot.csv")
colnames(ma)
ma <- ma[, c("clinical_diagnosis.x", "clinical_diagnosis2", "genetic_diagnosis", "solved_class", "gene_name")]
unique(ma$solved_class)
uniqueN(ma$genetic_diagnosis)


ma[ is.na(genetic_diagnosis) & solved_class == "Solved" , genetic_diagnosis := gene_name]
ma[ is.na(clinical_diagnosis.x)  , clinical_diagnosis.x := clinical_diagnosis2]
ma <- ma[, c("clinical_diagnosis.x",   "genetic_diagnosis" )]
ma <- ma[!duplicated(ma)]
ma <- ma[!is.na(genetic_diagnosis)]
ma <- ma[!is.na(clinical_diagnosis.x)]
unique(ma$clinical_diagnosis.x)
ma <- ma[ !(clinical_diagnosis.x %in% c("unspecified mitochondrial disorder", 
                                        "asymptomatic mutation carrier",
                                        "Unspecified mitochondrial disorder" , "Other" ))]

ma[genetic_diagnosis == "mtDNA mutation (unspecified)", genetic_diagnosis := "mtDNA mutation" ]
setnames(ma, c("clinical_diagnosis.x", "genetic_diagnosis"), c("clinical_diagnosis", "Gene_name" ))
ma <- ma[ !(clinical_diagnosis %in% c("Exclusion of mitochondrial disease" ))]

mito <- merge(mitog, ma , by  = "Gene_name", all = T)

iem <- fread("/Users/Mitya/Downloads/IEMbase_fz_old.xlsx - IEMbase.csv" )
setnames(iem, c("Gene Symbol"), c("Gene_name" ))

iem2 <- iem[Gene_name %in%  mitog$Gene_name]
nrow(iem2) / nrow(iem) * 100
miem <- merge(mito, iem , by  = "Gene_name", all.x = T)
miem <- miem[order(Gene_name, Ensembl_ID, IEMNosologyCode)]
write_csv(miem, "../mito_IEMbase.csv")

#####

df <- fread("/Users/Mitya/Downloads/mito_IEMbase - Sheet4.csv" ) 
df <- df[!duplicated(df)]

dfx <- merge(df, mitog,  by  = "Gene_name", all.x = T)
dfx$OMIM <- as.character(dfx$OMIM)
dfx[ is.na(OMIM), OMIM := OMIM_MIM_number]

dfx_miem <- merge(dfx[ , 1 : 5], iem , by  = "Gene_name", all.x = T)
dfx_miem <- dfx_miem[order(Syndrome, Gene_name, Ensembl_ID, IEMNosologyCode)]
write_csv(dfx_miem, "../selected_mito_IEMbase.csv")

uniqueN(iem$Gene_name)
uniqueN(dfx$Gene_name)
uniqueN(iem[ Gene_name %in% dfx$Gene_name]$Gene_name)

#############






###########################################
# Kremer Recon prot
###########################################
sag <- fread("/Users/Mitya/Downloads/4_RNAseq_annotation - sample_anno_transfered_gagneur_full.csv")
colnames(sag)
sag <- sag[ , c("FIBROBLAST_ID", "EXOME_ID", "RNA_ID", "KNOWN_MUTATION", "PROTEOME_ID" ,
                "BATCH" ,"SOLVED_BY",  "USE_IN_RNA_PAPER", "USE_COMMENT" , "OUTRIDER_GROUP", "COUNT_STRAND_SPECIFIC"  )]
sagK <- sag[(OUTRIDER_GROUP %like% "kremer") | is.na(RNA_ID)]
# sagKs <- sagK[is.na(PROTEOME_ID)]
# sagKs <- sagKs[ !is.na(KNOWN_MUTATION)]

# sagKs <- sagKs[!duplicated(FIBROBLAST_ID)]

prot <- fread("/Users/Mitya/Downloads/proteomics_annotation_QC.tsv" )


#sag[ , PID := paste0("P", FIBROBLAST_ID)]
#unique(sagKs$PID)
protS <- prot[Fibroblast_ID %in% sagK$FIBROBLAST_ID ]



protS <- protS[!is.na(Fibroblast_ID)]
protS <- protS[!(PROTEOME_ID %like% "Control_")]
protS <- protS[!(PROTEOME_ID %like% "Control")]
protS <- protS[!(PROTEOME_ID %like% "control")]
protS <- protS[ PLEX == "P11"]
protS <- protS[is.na(GROUP)]
sagX <- sagK[ FIBROBLAST_ID %in% protS$Fibroblast_ID ]


