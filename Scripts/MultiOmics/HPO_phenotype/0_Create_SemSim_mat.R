#'---
#' title: Similarity matrix
#' author: Dmitrii Smirnov
#' wb:
#'  input:
#'  - hpo: '`sm config["DATASETS"] + "/hp.obo"`'
#'  output:
#'  - tsm: '`sm config["DATASETS"] + "/hpGeneResnik.RDS"`'
#'  type: script
#'---

source("src/config.R")

# load HPO ontology
# hpo <- get_ontology("http://purl.obolibrary.org/obo/hp.obo", extract_tags="everything")

hpo <- get_ontology(snakemake@input$hpo)
# hpo <- get_ontology("datasets/hp.obo", extract_tags="everything")


# information content
ic <- descendants_IC(hpo)


#similarity of term pairs
tsm <- get_term_sim_mat(hpo, ic, method = "resnik")
#saveRDS(tsm, 'datasets/hpGeneResnik.RDS') 
saveRDS(tsm,  snakemake@output$tsm)


 







