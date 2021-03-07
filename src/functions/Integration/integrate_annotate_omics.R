##################################################
### OMICs integration and annotation functions ###
##################################################

require(data.table)
require(dplyr)
require(magrittr)




# Main integration function to combine rna seq and proteomics
combine_omics <- function(annotation = NULL,
                          RNA_ae = NULL, 
                          PROTEIN_ae = NULL, 
                          variants = NULL,
                          phenotype = NULL, 
                          variants_hpo = NULL,
                          Padj_threshold =  0.1,
                          Zscore_threshold = 3, 
                          l2FC_threshold = 1){
  
  ##########################################################
  # Combine RNA-seq and proteomics aberrant expression (AE)
  ##########################################################
  
  PROTEIN_ae[ , PROTEIN_INT := 2^PROTEIN_LOG2INT ]
  RNA_ae$geneID <- toupper(RNA_ae$geneID)
  rna_protein <-  merge(RNA_ae, PROTEIN_ae, by= c('SAMPLE_ID', 'geneID'), all= T) #outer join 
  rna_protein[, rank_rna := rank(normcounts +1, na.last= "keep"), by= geneID]
  rna_protein[, rank_protein := rank(PROTEIN_INT, na.last= "keep"), by= geneID]
  
  # Extended RNA and Protein outliar class
  rna_protein[, outlier_class := 'non_outlier']
  rna_protein[RNA_PADJ < Padj_threshold & PROTEIN_PADJ >= Padj_threshold, outlier_class := 'RNA-only']
  rna_protein[RNA_PADJ >= Padj_threshold & PROTEIN_PADJ < Padj_threshold, outlier_class := 'protein-only']
  
  rna_protein[RNA_PADJ < Padj_threshold & PROTEIN_PADJ < Padj_threshold, outlier_class := 'RNA-and-protein']
  
  rna_protein[RNA_PADJ < Padj_threshold & abs(PROTEIN_LOG2FC) > l2FC_threshold & PROTEIN_PVALUE < 0.05, outlier_class := 'RNA-and-protein']
  rna_protein[RNA_PADJ < Padj_threshold & abs(PROTEIN_ZSCORE) > Zscore_threshold , outlier_class := 'RNA-and-protein']
  
  rna_protein[ abs(RNA_LOG2FC) > l2FC_threshold & RNA_PVALUE < 0.05 & PROTEIN_PADJ < Padj_threshold, outlier_class := 'RNA-and-protein']
  rna_protein[ abs(RNA_ZSCORE) > Zscore_threshold & PROTEIN_PADJ < Padj_threshold, outlier_class := 'RNA-and-protein']
  
  
  # Add validation category
  rna_protein[, validated := F]
  rna_protein[outlier_class == 'Protein_outlier' , validated := T]
  rna_protein[PROTEIN_LOG2FC <= -1 & PROTEIN_PVALUE < 0.05 , validated := T]
  rna_protein[PROTEIN_ZSCORE <= -2 & PROTEIN_PVALUE < 0.05 , validated := T]
  
  
  # Add data availability flags
  rna_protein[ , RNA_seq_avaliable := SAMPLE_ID %in% unique(RNA_ae$SAMPLE_ID)]
  rna_protein[ , Proteomics_avaliable := SAMPLE_ID %in% unique(PROTEIN_ae$SAMPLE_ID)]
  
  # Add flag for gene products detected by RNA-seq or proteomics + passed the filters for AE
  rna_protein[, gene_detected := "RNA and protein detected" ]
  rna_protein[ is.na(RNA_ZSCORE), gene_detected := "no RNA" ]
  rna_protein[ is.na(PROTEIN_ZSCORE), gene_detected := "no Protein" ]
  rna_protein[ is.na(RNA_ZSCORE) & is.na(PROTEIN_ZSCORE), gene_detected := "no RNA and protein" ]
  
  
  
  ##########################################################
  # Add other OMICs data
  ##########################################################
  
  # Add variant annotation 
  if (!is.null(variants)) {
    rna_protein <- merge(variants, rna_protein, by= c('SAMPLE_ID', 'geneID'), all.y = T)
    
    # Add data availability flag
    rna_protein[ , WES_avaliable := SAMPLE_ID %in% unique(variants$SAMPLE_ID)]
    
    # Correct annotation
    rna_protein[WES_avaliable == T & is.na(rare), rare := F ]
    rna_protein[WES_avaliable == T & is.na(potential_biallelic), potential_biallelic := F ]
    
    # Add N variants info
    rna_protein[WES_avaliable == F , gene_class := "no data" ]
    rna_protein[WES_avaliable == T , gene_class := "no rare variant" ]
    rna_protein[WES_avaliable == T & rare == T , gene_class := "1 rare variant" ]
    rna_protein[WES_avaliable == T & potential_biallelic == T, gene_class := "rare pot. biallelic variants" ] 
    

    
    
    # UPD validation category - if variant information available, consider only genes with rare variants 
    rna_protein[WES_avaliable == T &  rare == F , validated := F]
    
  }
  
  # Add phenotype data    
  if (!is.null(phenotype)) {
    rna_protein <- merge(rna_protein, phenotype,  by= c('SAMPLE_ID', 'geneID'), all.x = T)
    
    # Add data availability flag
    rna_protein[ , HPO_avaliable := SAMPLE_ID %in% unique(phenotype$SAMPLE_ID)]
    
    # Mendelian disease gene flag
    rna_protein[ , MDG := geneID %in% unique(phenotype$geneID)]
    
    rna_protein[ is.na(Semantic_sim), Semantic_sim := 0]
  }
  
  
  # Add combined variants and phenotype data  
  if (!is.null(variants_hpo)) {
    rna_protein <- merge(variants_hpo, rna_protein, by= c('SAMPLE_ID', 'geneID'), all.y = T)
    
    # Add data availability flags
    rna_protein[ , WES_avaliable := SAMPLE_ID %in% unique(variants_hpo[WES_avaliable == T]$SAMPLE_ID)]
    rna_protein[ , HPO_avaliable := SAMPLE_ID %in% unique(variants_hpo[HPO_avaliable == T]$SAMPLE_ID)]
    
    # Mendelian disease gene flag
    rna_protein[ , MDG := geneID %in% unique(variants_hpo[!is.na(Semantic_sim)]$geneID)]
    rna_protein[ is.na(Semantic_sim), Semantic_sim := 0]
    
    # UPD annotation
    rna_protein[WES_avaliable == T & is.na(rare), rare := F ]
    rna_protein[WES_avaliable == T & is.na(potential_biallelic), potential_biallelic := F ]
    
    # UPD N variants info
    rna_protein[WES_avaliable == T & is.na( gene_class), gene_class := "no rare variant" ]
    rna_protein[WES_avaliable == F & is.na( gene_class), gene_class := "no data" ]
    
    # UPD validation category - if variant information available, consider only genes with rare variants 
    rna_protein[WES_avaliable == T &  rare == F , validated := F]
  }
  
  
  # Add causal gene annotation  
  if (!is.null(sa)) {
    diagnosed_cases <- sa[!is.na(KNOWN_MUTATION), c('SAMPLE_ID', 'KNOWN_MUTATION') ]
    diagnosed_cases[ , sample_gene := paste0(SAMPLE_ID, '_',  KNOWN_MUTATION )]
    
    rna_protein[ , sample_gene := paste0(SAMPLE_ID, '_',  geneID )]
    rna_protein[ , causal_gene := sample_gene %in% diagnosed_cases$sample_gene]
    
    # Add annotation for cases, where WES export was not possible or available
    rna_protein[causal_gene == T, rare := T ]
    
    rna_protein[ SAMPLE_ID == "OM25068" & geneID == "SLC25A4", rare := T ]
    rna_protein[ SAMPLE_ID == "OM25068" & geneID == "SLC25A4", potential_biallelic := T ]
    
    rna_protein[ SAMPLE_ID == "OM26649" & geneID == "AARS2", rare := T ]
    rna_protein[ SAMPLE_ID == "OM26649" & geneID == "AARS2", potential_biallelic := T ]
    
    rna_protein[ SAMPLE_ID == "OM65728" & geneID == "ACAD9", rare := T ]
    rna_protein[ SAMPLE_ID == "OM65728" & geneID == "ACAD9", potential_biallelic := T ]
    
    rna_protein[ SAMPLE_ID == "OM94976" & geneID == "DNAJC30", rare := T ]
    rna_protein[ SAMPLE_ID == "OM94976" & geneID == "DNAJC30", potential_biallelic := T ]
    
    # Correct annotation 
    rna_protein[rare == T , gene_class := "1 rare variant" ]
    rna_protein[ potential_biallelic == T, gene_class := "rare pot. biallelic variants" ] 
    
  }  
  return(rna_protein)
}



add_up_down_class <- function(omics,                           
                              Padj_threshold =  0.1,
                              Zscore_threshold = 3, 
                              l2FC_threshold = 1){
  
  #____Define over- under- expression outlier per omic__________________
  omics[, RNA_outlier_class := 'non_outlier']
  omics[RNA_PADJ < Padj_threshold & RNA_ZSCORE > 0 , RNA_outlier_class := 'RNA_overexpression']
  omics[RNA_PADJ < Padj_threshold & RNA_ZSCORE < 0 , RNA_outlier_class := 'RNA_underexpression']
  
  omics[, Protein_outlier_class := 'non_outlier']
  omics[ PROTEIN_PADJ < Padj_threshold & PROTEIN_ZSCORE > 0 , Protein_outlier_class := 'Protein_overexpression']
  omics[ PROTEIN_PADJ < Padj_threshold & PROTEIN_ZSCORE < 0 , Protein_outlier_class := 'Protein_underexpression']
  
  #____New RNA and Protein outlier over- under- expression class__________________
  
  omics[, up_down_outlier := 'non_outlier']
  
  omics[RNA_outlier_class == 'RNA_overexpression', up_down_outlier := 'RNA_overexpression']
  omics[RNA_outlier_class == 'RNA_underexpression', up_down_outlier := 'RNA_underexpression']
  
  omics[Protein_outlier_class == 'Protein_overexpression', up_down_outlier := 'Protein_overexpression']
  omics[Protein_outlier_class == 'Protein_underexpression', up_down_outlier := 'Protein_underexpression']
  
  omics[RNA_outlier_class == 'RNA_overexpression' & Protein_outlier_class == 'Protein_overexpression', up_down_outlier := 'RNA_Protein_overexpression']
  omics[RNA_outlier_class == 'RNA_underexpression' & Protein_outlier_class == 'Protein_underexpression', up_down_outlier := 'RNA_Protein_underexpression']
  
  omics[RNA_outlier_class == 'RNA_overexpression' &  PROTEIN_LOG2FC > l2FC_threshold & PROTEIN_PVALUE < 0.05, up_down_outlier := 'RNA_Protein_overexpression']
  omics[RNA_outlier_class == 'RNA_underexpression' & PROTEIN_LOG2FC < -l2FC_threshold & PROTEIN_PVALUE < 0.05, up_down_outlier := 'RNA_Protein_underexpression']
  
  
  omics[RNA_outlier_class == 'RNA_overexpression' &  PROTEIN_ZSCORE > Zscore_threshold, up_down_outlier := 'RNA_Protein_overexpression']
  omics[RNA_outlier_class == 'RNA_underexpression' & PROTEIN_ZSCORE < -Zscore_threshold, up_down_outlier := 'RNA_Protein_underexpression']
  
  omics[Protein_outlier_class == 'Protein_overexpression' & RNA_ZSCORE > Zscore_threshold , up_down_outlier := 'RNA_Protein_overexpression']
  omics[Protein_outlier_class == 'Protein_underexpression' & RNA_ZSCORE < -Zscore_threshold , up_down_outlier := 'RNA_Protein_underexpression']
  
  
  omics[Protein_outlier_class == 'Protein_overexpression' & 
          RNA_LOG2FC > l2FC_threshold & RNA_PVALUE < 0.05, up_down_outlier := 'RNA_Protein_overexpression']
  omics[Protein_outlier_class == 'Protein_underexpression' & 
          RNA_LOG2FC < -l2FC_threshold & RNA_PVALUE < 0.05, up_down_outlier := 'RNA_Protein_underexpression']
  
  
  
  omics[RNA_outlier_class == 'RNA_overexpression' & PROTEIN_ZSCORE < -Zscore_threshold, up_down_outlier := 'RNA_over_Protein_under']
  omics[Protein_outlier_class == 'Protein_underexpression' & RNA_ZSCORE > Zscore_threshold , up_down_outlier := 'RNA_over_Protein_under']
  
  omics[RNA_outlier_class == 'RNA_underexpression' & PROTEIN_ZSCORE > Zscore_threshold, up_down_outlier := 'RNA_under_Protein_over']
  omics[Protein_outlier_class == 'Protein_overexpression' & RNA_ZSCORE < -Zscore_threshold , up_down_outlier := 'RNA_under_Protein_over']
  
  
  omics[RNA_outlier_class == 'RNA_underexpression' & 
          Protein_outlier_class == 'Protein_overexpression', up_down_outlier := 'RNA_under_Protein_over']
  omics[RNA_outlier_class == 'RNA_underexpression' & 
          PROTEIN_LOG2FC > l2FC_threshold & PROTEIN_PVALUE < 0.05, up_down_outlier := 'RNA_under_Protein_over']
  omics[Protein_outlier_class == 'Protein_overexpression' & 
          RNA_LOG2FC < -l2FC_threshold & RNA_PVALUE < 0.05, up_down_outlier := 'RNA_under_Protein_over']
  
  
  omics[RNA_outlier_class == 'RNA_overexpression' & 
          Protein_outlier_class == 'Protein_underexpression', up_down_outlier := 'RNA_over_Protein_under']
  omics[RNA_outlier_class == 'RNA_overexpression' & 
          PROTEIN_LOG2FC < -l2FC_threshold & PROTEIN_PVALUE < 0.05, up_down_outlier := 'RNA_over_Protein_under']
  omics[Protein_outlier_class == 'Protein_underexpression' & 
          RNA_LOG2FC > l2FC_threshold & RNA_PVALUE < 0.05 , up_down_outlier := 'RNA_over_Protein_under']
  
  return( omics ) 
}


##########################################
### Add up- down- regulation outlier class
##########################################



add_up_down_class <- function(omics,                           
                              Padj_threshold =  0.1,
                              Zscore_threshold = 3, 
                              l2FC_threshold = 1){
  
  #____Define over- under- expression outlier per omic__________________
  omics[, RNA_outlier_class := 'non_outlier']
  omics[RNA_PADJ < Padj_threshold & RNA_ZSCORE > 0 , RNA_outlier_class := 'RNA_overexpression']
  omics[RNA_PADJ < Padj_threshold & RNA_ZSCORE < 0 , RNA_outlier_class := 'RNA_underexpression']
  
  omics[, Protein_outlier_class := 'non_outlier']
  omics[ PROTEIN_PADJ < Padj_threshold & PROTEIN_ZSCORE > 0 , Protein_outlier_class := 'Protein_overexpression']
  omics[ PROTEIN_PADJ < Padj_threshold & PROTEIN_ZSCORE < 0 , Protein_outlier_class := 'Protein_underexpression']
  
  #____New RNA and Protein outlier over- under- expression class__________________
  
  omics[, up_down_outlier := 'non_outlier']
  
  omics[RNA_outlier_class == 'RNA_overexpression', up_down_outlier := 'RNA_overexpression']
  omics[RNA_outlier_class == 'RNA_underexpression', up_down_outlier := 'RNA_underexpression']
  
  omics[Protein_outlier_class == 'Protein_overexpression', up_down_outlier := 'Protein_overexpression']
  omics[Protein_outlier_class == 'Protein_underexpression', up_down_outlier := 'Protein_underexpression']
  
  omics[RNA_outlier_class == 'RNA_overexpression' & Protein_outlier_class == 'Protein_overexpression', up_down_outlier := 'RNA_Protein_overexpression']
  omics[RNA_outlier_class == 'RNA_underexpression' & Protein_outlier_class == 'Protein_underexpression', up_down_outlier := 'RNA_Protein_underexpression']
  
  omics[RNA_outlier_class == 'RNA_overexpression' &  PROTEIN_LOG2FC > l2FC_threshold & PROTEIN_PVALUE < 0.05, up_down_outlier := 'RNA_Protein_overexpression']
  omics[RNA_outlier_class == 'RNA_underexpression' & PROTEIN_LOG2FC < -l2FC_threshold & PROTEIN_PVALUE < 0.05, up_down_outlier := 'RNA_Protein_underexpression']
  
  
  omics[RNA_outlier_class == 'RNA_overexpression' &  PROTEIN_ZSCORE > Zscore_threshold, up_down_outlier := 'RNA_Protein_overexpression']
  omics[RNA_outlier_class == 'RNA_underexpression' & PROTEIN_ZSCORE < -Zscore_threshold, up_down_outlier := 'RNA_Protein_underexpression']
  
  omics[Protein_outlier_class == 'Protein_overexpression' & RNA_ZSCORE > Zscore_threshold , up_down_outlier := 'RNA_Protein_overexpression']
  omics[Protein_outlier_class == 'Protein_underexpression' & RNA_ZSCORE < -Zscore_threshold , up_down_outlier := 'RNA_Protein_underexpression']
  
  
  omics[Protein_outlier_class == 'Protein_overexpression' & 
          RNA_LOG2FC > l2FC_threshold & RNA_PVALUE < 0.05, up_down_outlier := 'RNA_Protein_overexpression']
  omics[Protein_outlier_class == 'Protein_underexpression' & 
          RNA_LOG2FC < -l2FC_threshold & RNA_PVALUE < 0.05, up_down_outlier := 'RNA_Protein_underexpression']
  
  
  
  omics[RNA_outlier_class == 'RNA_overexpression' & PROTEIN_ZSCORE < -Zscore_threshold, up_down_outlier := 'RNA_over_Protein_under']
  omics[Protein_outlier_class == 'Protein_underexpression' & RNA_ZSCORE > Zscore_threshold , up_down_outlier := 'RNA_over_Protein_under']
  
  omics[RNA_outlier_class == 'RNA_underexpression' & PROTEIN_ZSCORE > Zscore_threshold, up_down_outlier := 'RNA_under_Protein_over']
  omics[Protein_outlier_class == 'Protein_overexpression' & RNA_ZSCORE < -Zscore_threshold , up_down_outlier := 'RNA_under_Protein_over']
  
  
  omics[RNA_outlier_class == 'RNA_underexpression' & 
          Protein_outlier_class == 'Protein_overexpression', up_down_outlier := 'RNA_under_Protein_over']
  omics[RNA_outlier_class == 'RNA_underexpression' & 
          PROTEIN_LOG2FC > l2FC_threshold & PROTEIN_PVALUE < 0.05, up_down_outlier := 'RNA_under_Protein_over']
  omics[Protein_outlier_class == 'Protein_overexpression' & 
          RNA_LOG2FC < -l2FC_threshold & RNA_PVALUE < 0.05, up_down_outlier := 'RNA_under_Protein_over']
  
  
  omics[RNA_outlier_class == 'RNA_overexpression' & 
          Protein_outlier_class == 'Protein_underexpression', up_down_outlier := 'RNA_over_Protein_under']
  omics[RNA_outlier_class == 'RNA_overexpression' & 
          PROTEIN_LOG2FC < -l2FC_threshold & PROTEIN_PVALUE < 0.05, up_down_outlier := 'RNA_over_Protein_under']
  omics[Protein_outlier_class == 'Protein_underexpression' & 
          RNA_LOG2FC > l2FC_threshold & RNA_PVALUE < 0.05 , up_down_outlier := 'RNA_over_Protein_under']
  
  return( omics ) 
}


##############
### Add OMIM columns
##############

add_omim_cols <- function(DT, gene_name_col = "gene_name", return_all_info = TRUE, pmim_link = T){
  omim_dt = readRDS("/s/project/mitoMultiOmics/db_data/omim-gene-pheno-cache.RDS")
  # omim_dt <- readRDS("../mitomultiomics/resource/omim_dt.Rds")
  omim_dt <- omim_dt[SYMBOL != "", .(SYMBOL, GMIM, PINH, PMIM)]
  omim_dt <- omim_dt[, SYMBOL := toupper(SYMBOL)]
  omim_dt <- omim_dt[, .SD[1], by = SYMBOL]   # In case of many GMIM, take the first only
  setnames(omim_dt, "GMIM", "OMIM")
  setnames(DT, gene_name_col, "gene_name")
  DT <- left_join(DT, omim_dt, by = c("gene_name" = "SYMBOL")) %>% as.data.table
  setnames(DT, "gene_name", gene_name_col)
  
  # Add a link to the OMIM website with the phenotype
  if(isTRUE(pmim_link)){
    DT[, OMIM_link := paste0('<a href=\"https://www.omim.org/entry/', PMIM,'\">', PMIM, '</a>'), by = 1:nrow(DT)]
    DT[is.na(PMIM), OMIM_link := NA]
    DT[grep(";", PMIM), OMIM_link := 'multiple_phenos']
  }
  
  
  # return_all_info: 
  ## TRUE: returns all the 3 columns
  ## FALSE: returns a T/F OMIM col
  if(isFALSE(return_all_info)){
    DT[, OMIM_gene := !is.na(PMIM) & PMIM != ""]
    DT[, c("OMIM", "PINH", "PMIM") := NULL]
  } 
  return(DT)
}

# gene_annot = add_omim_cols(gene_annot, gene_name_col = "gene_name")



add_nuclear_mito_DNA <- function(DT, gene_name_col = "gene_name"){
  ge <- readRDS("../genetic_diagnosis/resources/gencode.v19_with_gene_name.Rds")
  ge[, nDNA_mtDNA := "nDNA"]
  ge[seqnames == "MT", nDNA_mtDNA := "mtDNA"]
  setnames(DT, gene_name_col, "gene_name")
  DT <- left_join(DT, ge[,.(gene_name, nDNA_mtDNA)], by = "gene_name") %>% as.data.table
  setnames(DT, "gene_name", gene_name_col)
  return(DT)
}


add_gene_type <- function(DT, gene_name_col = "gene_name"){
  gt <- fread("/s/project/genetic_diagnosis/resource/gencode_v29_unique_gene_name.tsv")
  gt[, gene_name_unique := toupper(gene_name_unique)]
  setnames(DT, gene_name_col, "gene_name")
  DT <- left_join(DT, gt[,.(gene_name_unique, gene_type)], by = c("gene_name" = "gene_name_unique")) %>% as.data.table
  setnames(DT, "gene_name", gene_name_col)
  return(DT)
}




############################################
### add HGNC_symbol to UNIPROT ids
############################################

library(dplyr)

add_hgnc_symbol <- function(prot_table, hgnc_file, prot_id="UNIPROT_ID") {
  hgnc_mapping = fread(hgnc_file, header=TRUE)
  prot <- inner_join(hgnc_mapping[,c("UNIPROT_ID", "geneID" )], prot_table,  by = c("UNIPROT_ID" = prot_id)) %>% as.data.table()
  prot <- prot[prot$geneID !='',] # delete strange entries 
  prot <- prot[prot$UNIPROT_ID != '',  ]
  prot <- prot[!duplicated(prot ), ] # delete duplicated rows
  prot
}
