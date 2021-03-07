#'---
#' title: Patient organ system involvment
#' author: smirnovd
#' wb:
#'  input:
#'  - config: 'src/config.R'
#'  - sa: '`sm config["ANNOTATION"]`'
#'  - patient_hpo: '`sm config["RAW_HPO"]`'
#'  output:
#'  - patient_organs: '`sm config["PROC_DATA"] + "/HPO/Patients_affected_organs.tsv"`'
#'  - fig: '`sm config["FIGURE_DIR"] + "/Supplementary_figures/humans.pdf"`'
#' output: 
#'   html_document:
#'    code_folding: hide
#'    code_download: TRUE
#'---

source(snakemake@input$config)


# READ ANNOTATION
# sa <- fread('/s/project/mitoMultiOmics/multiOMICs_integration/raw_data/proteomics_annotation.tsv')
sa <- fread(snakemake@input$sa)
sa <- sa[USE_FOR_PROTEOMICS_PAPER == T]

# Load patient's HPO 
pat_hpo <- fread(snakemake@input$patient_hpo)
# pat_hpo <- fread('/s/project/mitoMultiOmics/multiOMICs_integration/raw_data/Patient_HPO_phenotypes.tsv')
pat_hpo <- pat_hpo[, c("SAMPLE_ID", "HPO_ID")]


# load HPO ontology
hpo <- get_ontology("http://purl.obolibrary.org/obo/hp.obo", extract_tags="everything")


# Extract HPO terms level 3 = organ systems
hpo_ID_to_hpo_term <- pat_hpo[, c("HPO_ID")] 
hpo_ID_to_hpo_term <- hpo_ID_to_hpo_term[!duplicated(hpo_ID_to_hpo_term)]
hpo_ID_to_hpo_term$umbrella_HPO_term <- lapply(hpo_ID_to_hpo_term$HPO_ID, function(hpo_id) {
  get_term_property(ontology=hpo, property="ancestors", term=hpo_id, as_names=TRUE)[3]
})
hpo_ID_to_hpo_term <- hpo_ID_to_hpo_term %>% 
  mutate_at('umbrella_HPO_term', paste)
#str(hpo_ID_to_hpo_term)
pat_hpo <- merge(pat_hpo, hpo_ID_to_hpo_term, by = "HPO_ID")


pat_hpo <- pat_hpo[, c("SAMPLE_ID", "umbrella_HPO_term")]
pat_hpo <-pat_hpo[!duplicated(pat_hpo ),]

sa <- sa[SAMPLE_ID %in% unique(pat_hpo$SAMPLE_ID)]
sa <- sa[!duplicated(sa )]

# create anatomy maps
library(gganatogram)                     
              
unique(pat_hpo$umbrella_HPO_term)
male <- as.data.table(hgMale_key) 
male$value <- NULL
male$sex <- rep('male', nrow(male))

female <- as.data.table(hgFemale_key)
female$value <- NULL
female$sex <- rep('female', nrow(female))
human <- rbind(male, female)     

human[type ==  "digestion", type:= "Abnormality of the digestive system" ]
human[type == "Abnormality of the digestive system" , colour := "orange"]

human[organ ==  "urinary_bladder", type:= "Abnormality of the genitourinary system" ]
human[type ==  "reproductive", type:= "Abnormality of the genitourinary system" ]
human[type == "Abnormality of the genitourinary system", colour := "yellow" ]

human[type ==  "nervous_system", type:= "Abnormality of the nervous system" ]
human[type == "Abnormality of the nervous system", colour := "lightblue" ]

human[organ ==  "skeletal_muscle", type:= "Abnormality of the musculoskeletal system"]
human[type == "Abnormality of the musculoskeletal system", colour := "red"]

# human[organ ==  "skeletal_muscle", type:= "Abnormality of the musculature"]
# human[organ ==  "bone" , type := "Abnormality of the skeletal system"  ]  
# human[type == "Abnormality of the skeletal system" , colour := "grey60" ]  


human[type ==  "circulation", type:= "Abnormality of the cardiovascular system"]
human[type == "Abnormality of the cardiovascular system", colour := "darkred"]

human[organ ==  "bone_marrow"  , type:= "Abnormality of blood and blood-forming tissues" ] 
human[type == "Abnormality of blood and blood-forming tissues" , colour := "darkorange"]  

human[type ==  "respiratory", type:= "Abnormality of the respiratory system"  ]
human[type == "Abnormality of the respiratory system" , colour := "steelblue" ]

human[organ ==  "lymph_node"   , type:= "Abnormality of the immune system"] 
human[type == "Abnormality of the immune system", colour := "darkgreen"] 

human[organ %in%  c("thyroid_gland", "pituitary_gland", "adrenal_gland")   , type:= "Abnormality of the endocrine system"]  
human[type == "Abnormality of the endocrine system" , colour := "green" ] 


# "Abnormality of the ear" 
# "Abnormality of the eye"
# "Abnormality of head or neck"
# "Growth abnormality"  
# "Abnormality of limbs" 
# "Abnormality of metabolism/homeostasis"
# "Onset"  
# "Abnormality of the integument" 
# "Triggered by"

human <- human[type !=  "other"]
human <- human[organ !=  "tongue"]
human <- human[organ !=  "throat"]
human <- human[organ !=  "trachea"]
human <- human[organ !=  "esophagus"]
human <- human[organ !=  "vas_deferens"]
human <- human[organ !=  "nerve"] # too many 
human <- human[organ !=  "spinal_cord"] # too many 



patients <- merge(sa[ , c("SAMPLE_ID", "gender") ], pat_hpo, by= "SAMPLE_ID" )
setnames(patients, c("gender", "umbrella_HPO_term"), c("sex", "type"))

patients <- merge(patients, human, by= c( "sex", "type"), allow.cartesian=TRUE)


# Save results
# write_tsv(patients, '/s/project/mitoMultiOmics/multiOMICs_integration/processed_data/HPO/Patients_affected_organs.tsv')
write_tsv(patients,  snakemake@output$patient_organs)


p_male <- gganatogram(data=  patients[sex == "male"], fillOutline='white', organism='human', sex= "male", fill="colour") + theme_void()
p_female <- gganatogram(data=  patients[sex == "female"], fillOutline='white', organism='human', sex= "female", fill="colour")+ theme_void()
humans <- p_male + p_female 

#+ fig.width=10, fig.height=7
humans


pdf(snakemake@output$fig, # "/s/project/mitoMultiOmics/multiOMICs_integration/Figures/humans.pdf",  
    width = 10, height = 7,  useDingbats=FALSE )
print(humans) 
dev.off()
