#'---
#' title: Gene semantic similarity
#' author: smirnovd
#' wb:
#'  input:
#'  - sa: '`sm config["ANNOTATION"]`'
#'  - patient_hpo: '`sm config["RAW_HPO"]`'
#'  output:
#'  - semantic_similariy: '`sm config["PROC_DATA"] + "/HPO/Patient_Gene_semantic_similariy.tsv"`'
#'  threads: 40
#' output:
#'   html_document:
#'    code_folding: hide
#'    code_download: TRUE
#'---

source("src/config.R")


# READ ANNOTATION
sa <- fread(snakemake@input$sa)
sa <- sa[USE_FOR_PROTEOMICS_PAPER == T]


# load HPO ontology
# data(hpo) 
hpo <- get_ontology("http://purl.obolibrary.org/obo/hp.obo", extract_tags="everything")


# Load patient's HPO 
pat_hpo <- fread(snakemake@input$patient_hpo)


pat_hpo <- pat_hpo[, c("SAMPLE_ID", "HPO_ID")]
# renmove duplicates from HPOs per patient 
pat_hpo <-pat_hpo[!duplicated(pat_hpo ),]
pat_hpo <-pat_hpo[HPO_ID %in%  unique(hpo$id)]


hpo_gene <- fread("http://compbio.charite.de/jenkins/job/hpo.annotations/lastSuccessfulBuild/artifact/util/annotation/phenotype_to_genes.txt")
colnames(hpo_gene) <- c("HPO_ID", "HPO_name", "entrezID", "geneID", "Additional_Info",  "source", "disease-ID")
hpo_gene <- hpo_gene[, c("geneID", "HPO_ID")]
hpo_gene <- hpo_gene[!duplicated(hpo_gene)]
hpo_gene <-hpo_gene[HPO_ID %in%  unique(hpo$id)]


# Fix gene names
hpo_gene[geneID == "YARS1", geneID := "YARS" ]
hpo_gene[geneID == "MICOS13", geneID := "C19ORF70" ]
hpo_gene[geneID == "ADPRS", geneID := "ADPRHL2" ]


# Add known genotype - phenotype associations 
hpo_gene <- rbind(hpo_gene, data.table( geneID = "DNAJC30", HPO_ID = c("HP:0002076", "HP:0000648", "HP:0000505", "HP:0011462", "HP:0001112", "HP:0003581", "HP:0001263", "HP:0001332", "HP:0001251", "HP:0002134", "HP:0011923", "HP:0011463"  )  ) )
hpo_gene <- rbind(hpo_gene, data.table( geneID = "MRPS25", HPO_ID = c("HP:0003581", "HP:0011445", "HP:0001274", "HP:0002421", "HP:0001266", "HP:0002509", "HP:0001348", "HP:0001274", "HP:0007333", "HP:0001324", 
                                                                      "HP:0001385", "HP:0000846", "HP:0004322", "HP:0000252", "HP:0011924", "HP:0008347", "HP:0001270", "HP:0003676", "HP:0012378", "HP:0002015")  ) ) 
hpo_gene <- rbind(hpo_gene, data.table( geneID = "TXNIP", HPO_ID = c("HP:0004902", "HP:0001943", "HP:0003355", "HP:0006568", "HP:0003658", "HP:0004322", "HP:0001508","HP:0001252", "HP:0011968")  ) )
hpo_gene <- rbind(hpo_gene, data.table( geneID = "YARS", HPO_ID = c("HP:0001510", "HP:0001263", "HP:0007305", "HP:0000407", "HP:0000639", "HP:0002611", "HP:0001738", "HP:0001943", "HP:0001903", "HP:0000093", "HP:0100806", "HP:0006528")  ) )

hpo_gene <- hpo_gene[!duplicated(hpo_gene)]
sa <- sa[SAMPLE_ID %in% unique(pat_hpo$SAMPLE_ID) ]

# extract ancestors
hp_ancestors <- hpo$ancestors
# create HPO named list of genes
info <- unstack(hpo_gene, geneID~HPO_ID)
hpByGene <- unstack(hpo_gene, HPO_ID~geneID)

# information content
ic <- descendants_IC(hpo)


## Compute semantic similarity between HP of interest and all HP terms
## This step is time consumming and can be parallelized.
SSscores <- data.frame()
for (patient in unique(sa$SAMPLE_ID)){
  hpOfInterest <- pat_hpo[SAMPLE_ID == patient]$HPO_ID
  
  hpGeneResnik <- compareHPSets(
    hpSet1=names(ic), hpSet2=hpOfInterest,
    IC=ic,
    ancestors=hp_ancestors,
    method="Resnik",
    #BPPARAM= MulticoreParam(40)
    BPPARAM= MulticoreParam(snakemake@threads)
  )
  
  ## Group the results by gene
  hpMatByGene <- lapply(
    hpByGene,
    function(x){
      hpGeneResnik[x, , drop=FALSE]
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
  res <- as.data.table(res[, c("SAMPLE_ID", "geneID", "Semantic_sim" )])
  res <- res[order(Semantic_sim, decreasing = T)]
  res$Rank_SSs <- seq(1, nrow(res))
  SSscores <- rbind(SSscores, res)
}

hist(SSscores$Semantic_sim, xlab = "Semantic similariy score",   
      main = "Histogram of Semantic similariy scores", col = "gray80", breaks = 30)

write_tsv(SSscores,  snakemake@output$semantic_similariy)

