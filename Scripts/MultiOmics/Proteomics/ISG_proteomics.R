sa <- fread('kaisa_prep/raw_data/proteomics_annotation.tsv')
uniqueN(sa$PROTEOME_ID)
sa[ , .N, by = GROUP]
sa <- sa[FAILED == F]
saU <- sa[USE == T]
protein_annotation <- fread("kaisa_prep/raw_data/gene_annotation_v29.tsv")

prot <- readRDS(  'kaisa_prep/processed_data/limma_results.rds')

rna <- readRDS(  '/Users/Mitya/Downloads/OUTRIDER_results_all_fib_ss.Rds')
head(rna)
uniqueN(rna$sampleID)

pa <- protein_annotation[ gene_name %in% c( "IFI44L", "IFI27", "USP18", "IFI6")]
rISG <- rna[ geneID %in% pa$gene_id_unique]

rISG[ , ISG := mean( zScore , na.rm = T), by = sampleID]
rISG <- rISG[ , c("sampleID", "ISG" )]
rISG <- rISG[ !duplicated( rISG)]
rISGx <- merge(rISG, sad[ , c("RNA_ID" , "Sex", "DISEASE", "pdescription")] , 
               by.x = "sampleID", by.y = "RNA_ID")
sad <- fread('/Users/Mitya/Downloads/drop_sa.tsv')
ggplot( rISG)

ggplot(rISGx[ Sex %in% c("male", "female")], aes(Sex, ISG)) +
  geom_quasirandom(aes( colour = Sex)) +
  stat_compare_means() +
  scale_color_ptol()+
  theme_classic()

