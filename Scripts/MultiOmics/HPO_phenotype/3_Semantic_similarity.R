#'---
#' title: Phjenotype semantic similarity
#' author: smirnovd
#' wb:
#'  input:
#'  - sa: '`sm config["ANNOTATION"]`'
#'  - patient_hpo: '`sm config["RAW_HPO"]`'
#'  - hpo_gene: '`sm config["RAW_DATA"] + "/hpo_gene.txt"`'
#'  - tsm: '`sm config["DATASETS"] + "/hpGeneResnik.RDS"`'
#'  output:
#'  - semantic_similariy: '`sm config["PROC_DATA"] + "/HPO/Patient_Gene_semantic_similariy.tsv"`'
#'  threads: 4
#' output:
#'   html_document:
#'    code_folding: hide
#'    code_download: TRUE
#'---

source("src/config.R")

threads <- snakemake@threads

# READ ANNOTATION
sa <- fread(snakemake@input$sa)
sa <- sa[USE_FOR_PROTEOMICS_PAPER == T]


# load HPO ontology
# hpo <- get_ontology("http://purl.obolibrary.org/obo/hp.obo", extract_tags="everything")
hpo <- get_ontology("datasets/hp.obo", extract_tags="everything")

# load precomputed matrix of pairwise HPO term similarities
#hpGeneResnik <- readRDS('datasets/hpGeneResnik.RDS' )
hpGeneResnik <- readRDS(snakemake@input$tsm )



# Load patient's HPO 
pat_hpo <- fread(snakemake@input$patient_hpo)
pat_hpo <- pat_hpo[, c("SAMPLE_ID", "HPO_ID")]
pat_hpo <-pat_hpo[!duplicated(pat_hpo ),]

sa <- sa[SAMPLE_ID %in% unique(pat_hpo$SAMPLE_ID) ]


#load gene - hpo 
hpo_gene <- fread(snakemake@input$hpo_gene)
# hpo_gene <- fread("../omicsDagnostics_data/raw_data/hpo_gene.txt")

gene_anno <- hpo_gene[ , c("geneID",  "type")]
gene_anno <- gene_anno[!duplicated(gene_anno)]


hpo_gene <- hpo_gene[, c("geneID", "HPO_ID")]
hpo_gene <- hpo_gene[!duplicated(hpo_gene)]




  


# create HPO named list of genes
info <- unstack(hpo_gene, geneID~HPO_ID)
hpByGene <- unstack(hpo_gene, HPO_ID~geneID)

# information content
ic <- descendants_IC(hpo)


## Compute semantic similarity between HP of interest and all HP terms
## This step is time consuming and can be paralleled.

hpGeneResnik <- hpGeneResnik[ , unique(pat_hpo$HPO_ID)]


# patient <-  "103527"
SSscores <- data.frame()
for ( patient in unique(sa$SAMPLE_ID) )  { #    %dopar%   patient in unique(pat_hpo$PHENOME_ID)
  
  # start_time <- Sys.time()
  hpOfInterest <- pat_hpo[SAMPLE_ID == patient]$HPO_ID
  
 
  # Extract similarity matrix per patient 
  hpGeneResnik_Pat <- hpGeneResnik[ , hpOfInterest]
  
  
  ## Group the results by gene
  hpMatByGene <- lapply(
    hpByGene,  
    function(x){
      hpGeneResnik_Pat[x, , drop=FALSE]
    }
  )
  
  ## Compute the corresponding scores
  Semantic_sim <- unlist(lapply(
    hpMatByGene,
    hpSetCompSummary,
    method="bma", direction="symSim"
  ))
  
  res <- as.data.frame(Semantic_sim)
  res$geneID <- rownames(res) 
  res$SAMPLE_ID <- rep(  patient, nrow(res))
  
  
  sum_ic <- sum(ic[ get_ancestors(hpo, hpOfInterest) ]) # hpOfInterest 
  res$sumIC <- rep( sum_ic , nrow(res)) 
  
  res <- as.data.table(res)
  res <- res[order(Semantic_sim, decreasing = T)]
  
  # Add re-scaling
  res[ , SemSim := scales:::rescale(Semantic_sim, to = c(1, 100))]
  res <- res[order(SemSim, decreasing = T)]  
  res$Rank_SemSim <- seq(1, nrow(res))
  
  # Add rankings per type
  resA <- merge(res, gene_anno, by = "geneID", all.x = T)
  resA <- resA[!duplicated(resA)]
  
  
  res_gene <- resA[ type == "gene"]
  res_gene <- res_gene[order(SemSim, decreasing = T)]
  res_gene$Rank_SemSim_type <- seq(1, nrow(res_gene))
  
  res_pat <- resA[ type == "patient"]
  res_pat <- res_pat[order(SemSim, decreasing = T)]  
  res_pat$Rank_SemSim_type <- seq(1, nrow(res_pat))
  
  
  res <- rbind(res_gene, res_pat)
  res <- res[!duplicated(res)]
  res <- res[order(SemSim, decreasing = T)]
  
  # stop_time <- Sys.time()
  # ana_time <- stop_time - start_time
  # 
  
  SSscores <- rbind(SSscores, res)
  # cat(patient, "done in", ana_time, "\n")
  
}


hist(SSscores$SemSim, xlab = "Semantic similariy score",   
      main = "Histogram of Semantic similariy scores", col = "gray80", breaks = 30)

write_tsv(SSscores,  snakemake@output$semantic_similariy)

