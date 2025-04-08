# Load required libraries
library(data.table)
library(yaml)

# Load configuration
config <- yaml::yaml.load_file('wbuild.yaml')
ANNOTATION <- config$ANNOTATION
RAW_HPO <- config$RAW_HPO
PROC_DATA <- config$PROC_DATA

# Process sample annotation
sa <- fread(ANNOTATION)
sa <- sa[USE_FOR_PROTEOMICS_PAPER == TRUE]
sa[, solved_class := "unsolved"]
sa[CATEGORY %in% c("IIb", "IIc"), solved_class := "candidate"]
sa[CATEGORY %in% c("I.m", "I", "IIa", "III"), solved_class := "solved"]

# Remove unsolved cases
sa <- sa[solved_class != "unsolved"]

# Load integrated omics data
rp <- readRDS(paste0(PROC_DATA, "/integration/patient_omics_full.RDS")) %>% as.data.table()

# Handle NA values in semantic similarity
rp[is.na(Semantic_sim), Semantic_sim := 0]

# Rescale semantic similarity to 0-100 per sample
rp[, Semantic_sim_scaled := 0]
rp[Semantic_sim > 0, Semantic_sim_scaled := (Semantic_sim - min(Semantic_sim)) / (max(Semantic_sim) - min(Semantic_sim)) * 100, by = SAMPLE_ID]

# Clean and prepare data
rp[, sample_gene := paste0(SAMPLE_ID, "_", geneID)]

# Subset samples in the annotation
sa <- sa[SAMPLE_ID %in% unique(rp$SAMPLE_ID)]


rp <- rp[SAMPLE_ID %in% unique(sa$SAMPLE_ID) ]

# Load complexes data
complexes <- readRDS(paste0(PROC_DATA, "/Complexes/Complex_outliers_LIMMA.rds")) %>% as.data.table()
complexes <- complexes[SAMPLE_ID %in% sa$SAMPLE_ID]

# Load and process patient HPO data
pat_hpo <- fread(RAW_HPO)
pat_hpo <- pat_hpo[SAMPLE_ID %in% sa$SAMPLE_ID]


hpo <- get_ontology("datasets/hp.obo", extract_tags="everything")
hpo_terms <- data.table(
  HPO_ID = hpo$id,
  HPO_term = hpo$name
)

pat_hpo <- merge(pat_hpo,hpo_terms, by = "HPO_ID" )


# Load patients affected organs
patients <- fread(paste0(PROC_DATA, '/HPO/Patients_affected_organs.tsv'))
patients <- patients[SAMPLE_ID %in% sa$SAMPLE_ID]

# Create shiny_data directory if it doesn't exist
dir.create("shiny_data", showWarnings = FALSE)

# Save processed data
saveRDS(sa, "omicsDiagnosticsAPP/shiny_data/sample_annotation.RDS")
saveRDS(rp, "omicsDiagnosticsAPP/shiny_data/patient_omics.RDS")
saveRDS(pat_hpo, "omicsDiagnosticsAPP/shiny_data/patient_hpo.RDS")
saveRDS(patients, "omicsDiagnosticsAPP/shiny_data/patients_organs.RDS")
saveRDS(complexes, "omicsDiagnosticsAPP/shiny_data/complexes.RDS")
 
# Create minimal config file
minimal_config <- list(
  ANNOTATION = "shiny_data/sample_annotation.RDS",
  PROC_DATA = "shiny_data"
)

# Save minimal config
yaml::write_yaml(minimal_config, "shiny_data/config.yaml") 