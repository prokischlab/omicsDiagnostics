
#load functions
source("kaisa_prep/Scripts/config.R")

total_protein_mass <- 15000 # 15 µg , 15000 is in nanograms (ng)

# Function to process MaxQuant proteinGroups.txt data, with options for original or absolute protein quantification
process_proteinGroups <- function(data, 
                                  approach = c("ProteinIDs", "MajorityProteinIDs"), 
                                  total_protein_mass = 15000, 
                                  output_type = c("intensity", "absolute")) {
  
  # Ensure approach and output_type are valid
  approach <- match.arg(approach)
  output_type <- match.arg(output_type)
  
  # Convert input data to data.table for better performance
  setDT(data)
  
  # Step 1: Filter out contaminants and reverse hits
  filtered_data <- data[!(Reverse == "+" | `Potential contaminant` == "+")]
  
  # Step 2: Extract relevant Ensembl Gene IDs (gene_id_unique) for aggregation
  if (approach == "ProteinIDs") {
    # Split 'Protein IDs' into individual rows
    processed_data <- filtered_data[, .(Protein_IDs = unlist(strsplit(`Protein IDs`, ";"))), by = "Protein IDs"]
    
    # Extract ENSG (gene) ID for aggregation, including the version and underscore part
    processed_data[, gene_id_unique := str_extract(Protein_IDs, "ENSG[0-9]+(?:\\.[0-9]+)?(?:_[0-9]+)?")]
    
    # Merge processed_data back with filtered_data based on `Protein IDs` to keep the intensity reporter columns
    processed_data <- merge(processed_data, filtered_data, by = "Protein IDs", all.x = TRUE)
    
  } else if (approach == "MajorityProteinIDs") {
    # Split 'Majority protein IDs' into individual rows
    processed_data <- filtered_data[, .(Majority_Protein_IDs = unlist(strsplit(`Majority protein IDs`, ";"))), by = "Protein IDs"]
    
    # Extract ENSG (gene) ID for aggregation, including the version and underscore part
    processed_data[, gene_id_unique := str_extract(Majority_Protein_IDs, "ENSG[0-9]+(?:\\.[0-9]+)?(?:_[0-9]+)?")]
    
    # Merge processed_data back with filtered_data based on `Protein IDs` to keep the intensity reporter columns
    processed_data <- merge(processed_data, filtered_data, by = "Protein IDs", all.x = TRUE)
  }
  
  # Step 3: If the output is "absolute", check for protein size adjustment factors
  if (output_type == "absolute") {
    if ("iBAQ peptides" %in% colnames(processed_data)) {
      processed_data[, Adjusted_Intensity_Factor := `iBAQ peptides`]
    } else if ("Sequence length" %in% colnames(processed_data)) {
      processed_data[, Adjusted_Intensity_Factor := `Sequence length`]
    } else {
      # If neither iBAQ peptides nor sequence length is available, return intensity intensities by default
      warning("No 'iBAQ peptides' or 'Sequence length' column available for adjustment. Returning intensity intensities.")
      return(processed_data[, c("gene_id_unique", grep("Reporter intensity corrected", colnames(processed_data), value = TRUE)), with = FALSE])
    }
    
    processed_data[ Adjusted_Intensity_Factor == 0 , Adjusted_Intensity_Factor := NA]
    
    # Step 4: Adjust reporter intensities for each protein based on protein size
    reporter_columns <- grep("Reporter intensity corrected", colnames(processed_data), value = TRUE)
    
    for (col in reporter_columns) {
      processed_data[, (paste0("Adjusted_", col)) := get(col) / Adjusted_Intensity_Factor]
    }
    
    # Step 5: Aggregate proteins at the gene (gene_id_unique) level
    aggregated_data <- processed_data[, lapply(.SD, sum, na.rm = TRUE), by = gene_id_unique, .SDcols = paste0("Adjusted_", reporter_columns)]
    
    # Step 6: Scale adjusted intensities to absolute protein quantities using the total protein mass
    for (i in seq_along(reporter_columns)) {
      adjusted_col <- paste0("Adjusted_", reporter_columns[i])
      original_col <- reporter_columns[i]
      aggregated_data[, (original_col) := (get(adjusted_col) / sum(get(adjusted_col), na.rm = TRUE)) * total_protein_mass]
    }
    
    # Step 7: Return the gene_id_unique and the columns with absolute protein quantifications
    return(aggregated_data[, c("gene_id_unique", reporter_columns), with = FALSE])
  }
  
  # If output_type is "intensity", return the original intensities
  return(processed_data[, c("gene_id_unique", grep("Reporter intensity corrected", colnames(processed_data), value = TRUE)), with = FALSE])
}


# min(processed_data$Adjusted_Intensity_Factor)
# 
# b26 <- paste0( "Adjusted_", sa[Experiment == "B26_mix1" ]$ID_raw_new)
# View(processed_data[ , ..b26])
# df <- aggregated_data[ , ..b26]
# sum(df$`Adjusted_Reporter intensity corrected 1 B26_mix1`, na.rm = T)
# aggregated_data[ , ..b26]
# > sum(aggregated_data$`Adjusted_Reporter intensity corrected 1 B26_mix1`, na.rm = T)
# [1] Inf

gencode29_prot <- fread("kaisa_prep/raw_data/gene_annotation_v29.tsv")

sa <- fread('kaisa_prep/raw_data/proteomics_annotation.tsv')
sa[ Experiment == "B37_Mix6", Experiment := "B_37_Mix6" ] 
sa[ , ID_raw_new := paste("Reporter intensity corrected", Fraction, Experiment)]

#########################################################
############# Ensemble proteomics #######################
#########################################################

prot_V29 <-  fread("kaisa_prep/raw_data/combined_2023_03_24_TMT_11plex_proteinGroups.txt")
#reporter_columns <- grep("Reporter intensity corrected", colnames(prot_V29), value = TRUE)
raw_prot <- process_proteinGroups(prot_V29, approach = "ProteinIDs", total_protein_mass = 15000, output_type = "absolute")
b26 <- paste0(  sa[Experiment == "B26_mix1" ]$ID_raw_new)
View(raw_prot[ , ..b26])
setnames(raw_prot, sa$ID_raw_new, sa$PROTEOME_ID, skip_absent=TRUE)
rm(prot_V29)
raw_prot <- raw_prot[!is.na(gene_id_unique)]

prot_new <-  fread("kaisa_prep/raw_data/20231110_P014_34_proteinGroups.txt")
prot_newA <- process_proteinGroups(prot_new, approach = "ProteinIDs", total_protein_mass = 15000, output_type = "absolute")
setnames(prot_newA, sa$ID_raw_new, sa$PROTEOME_ID, skip_absent=TRUE)
prot_newA <- prot_newA[!is.na(gene_id_unique)]
rm(prot_new)
np <- (colnames(prot_newA)[!(colnames(prot_newA) %in% colnames(raw_prot)[-1]) ])
raw_prot <- merge(raw_prot, prot_newA[ ,  ..np ], by = "gene_id_unique", all = T)
rm(prot_newA)
raw_prot <- raw_prot[!is.na(gene_id_unique)]

prot_new <-  fread("kaisa_prep/raw_data/20231218_P014_35_proteinGroups.txt")
prot_newA <- process_proteinGroups(prot_new, approach = "ProteinIDs", total_protein_mass = 15000, output_type = "absolute")
setnames(prot_newA, sa$ID_raw_new, sa$PROTEOME_ID, skip_absent=TRUE)
prot_newA <- prot_newA[!is.na(gene_id_unique)]
rm(prot_new)
np <- (colnames(prot_newA)[!(colnames(prot_newA) %in% colnames(raw_prot)[-1]) ])
raw_prot <- merge(raw_prot, prot_newA[ ,  ..np ], by = "gene_id_unique", all = T)
rm(prot_newA)
raw_prot <- raw_prot[!is.na(gene_id_unique)]

prot_new <-  fread("kaisa_prep/raw_data/20231218_P014_36_proteinGroups.txt")
prot_newA <- process_proteinGroups(prot_new, approach = "ProteinIDs", total_protein_mass = 15000, output_type = "absolute")
setnames(prot_newA, sa$ID_raw_new, sa$PROTEOME_ID, skip_absent=TRUE)
prot_newA <- prot_newA[!is.na(gene_id_unique)]
rm(prot_new)
np <- (colnames(prot_newA)[!(colnames(prot_newA) %in% colnames(raw_prot)[-1]) ])
raw_prot <- merge(raw_prot, prot_newA[ ,  ..np ], by = "gene_id_unique", all = T)
rm(prot_newA)
raw_prot <- raw_prot[!is.na(gene_id_unique)]

prot_new <-  fread("kaisa_prep/raw_data/20240209_P014_37m7_proteinGroups.txt")
prot_newA <- process_proteinGroups(prot_new, approach = "ProteinIDs", total_protein_mass = 15000, output_type = "absolute")
setnames(prot_newA, sa$ID_raw_new, sa$PROTEOME_ID, skip_absent=TRUE)
prot_newA <- prot_newA[!is.na(gene_id_unique)]
rm(prot_new)
np <- (colnames(prot_newA)[!(colnames(prot_newA) %in% colnames(raw_prot)[-1]) ])
raw_prot <- merge(raw_prot, prot_newA[ ,  ..np ], by = "gene_id_unique", all = T)
rm(prot_newA)
raw_prot <- raw_prot[!is.na(gene_id_unique)]

prot_new <-  fread("kaisa_prep/raw_data/20240307_P014_37_m3_4_5_proteinGroups.txt")
prot_newA <- process_proteinGroups(prot_new, approach = "ProteinIDs", total_protein_mass = 15000, output_type = "absolute")
setnames(prot_newA, sa$ID_raw_new, sa$PROTEOME_ID, skip_absent=TRUE)
prot_newA <- prot_newA[!is.na(gene_id_unique)]
rm(prot_new)
np <- (colnames(prot_newA)[!(colnames(prot_newA) %in% colnames(raw_prot)[-1]) ])
raw_prot <- merge(raw_prot, prot_newA[ ,  ..np ], by = "gene_id_unique", all = T)
rm(prot_newA)
raw_prot <- raw_prot[!is.na(gene_id_unique)]

prot_new <-  fread("kaisa_prep/raw_data/20240308_P014_37_m1_2_6_8_proteinGroups.txt")
prot_newA <- process_proteinGroups(prot_new, approach = "ProteinIDs", total_protein_mass = 15000, output_type = "absolute")
setnames(prot_newA, sa$ID_raw_new, sa$PROTEOME_ID, skip_absent=TRUE)
prot_newA <- prot_newA[!is.na(gene_id_unique)]
rm(prot_new)
np <- (colnames(prot_newA)[!(colnames(prot_newA) %in% colnames(raw_prot)[-1]) ])
raw_prot <- merge(raw_prot, prot_newA[ ,  ..np ], by = "gene_id_unique", all = T)
rm(prot_newA)
raw_prot <- raw_prot[!is.na(gene_id_unique)]

prot_new <-  fread("kaisa_prep/raw_data/20240808_P014_38_m1_proteinGroups.txt")
prot_newA <- process_proteinGroups(prot_new, approach = "ProteinIDs", total_protein_mass = 15000, output_type = "absolute")
setnames(prot_newA, sa$ID_raw_new, sa$PROTEOME_ID, skip_absent=TRUE)
prot_newA <- prot_newA[!is.na(gene_id_unique)]
rm(prot_new)
np <- (colnames(prot_newA)[!(colnames(prot_newA) %in% colnames(raw_prot)[-1]) ])
raw_prot <- merge(raw_prot, prot_newA[ ,  ..np ], by = "gene_id_unique", all = T)
rm(prot_newA)
raw_prot <- raw_prot[!is.na(gene_id_unique)]
raw_prot[ is.na(raw_prot)] <- 0

write_tsv(raw_prot, "kaisa_prep/raw_data/proteomics_not_normalized_absolute_ens.tsv")
raw_prot <- fread("kaisa_prep/raw_data/proteomics_not_normalized_absolute_ens.tsv")

# Aggregate on gene level ####################################################################################################
prot_new_a <- merge(gencode29_prot[, c("gene_id_unique", "gene_name")],  raw_prot, by = "gene_id_unique" )
prot_new_a$gene_id_unique <- NULL
prot_new_a[ duplicated(gene_name), ]$gene_name
setnames(prot_new_a, "gene_name", "geneID")

# Aggregate proteins at the gene (gene_id_unique) level
prot_new_agr <- prot_new_a[, lapply(.SD, sum, na.rm = TRUE), by = geneID ]
# View( prot_new_a[geneID == "SOD2" ] )


prot_new <- as.data.frame(prot_new_agr)
# rownames(prot_new) <- prot_new$geneID
# prot_new$geneID <- NULL


write_tsv(prot_new_agr, "kaisa_prep/raw_data/proteomics_not_normalized_absolute.tsv")



