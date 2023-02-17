#'---
#' title: HPO match Fresard et al 2019  
#' author: smirnovd
#' wb:
#'  input:
#'  - config: 'src/config.R'
#'  - sa: '`sm config["ANNOTATION"]`'
#'  - patient_hpo: '`sm config["RAW_HPO"]`'
#'  output:
#'  - hpo_match: '`sm config["PROC_DATA"] + "/HPO/Patients_HPO_Gene_mapping.tsv"`'
#' output:
#'   html_document:
#'    code_folding: hide
#'    code_download: TRUE
#'---


source(snakemake@input$config)


# READ ANNOTATION
sa <- fread(snakemake@input$sa)
sa <- sa[USE_FOR_PROTEOMICS_PAPER == T]

# Load patient's HPO 
pat_hpo <- fread(snakemake@input$patient_hpo)


pat_hpo <- pat_hpo[, c("SAMPLE_ID", "HPO_ID")]
pat_hpo <-pat_hpo[!duplicated(pat_hpo )]




# HPO mapping
# load HPO ontology
hpo <- get_ontology("http://purl.obolibrary.org/obo/hp.obo", extract_tags="everything")




# extract parent and chiled for a given term 
hpo_names <- data.frame(HPO_ID = hpo$id, HPO_name = hpo$name,  stringsAsFactors = F) #save names

# extract parents
parents <- ldply(hpo$parents, data.frame, stringsAsFactors = F)
colnames(parents) <- c("HPO_ID", "relative")
parents$relation <- "parent"

# extract children
children <- ldply(hpo$children, data.frame, stringsAsFactors = F)
colnames(children) <- c("HPO_ID", "relative")
children$relation <- "child"

# combine parents and children
hpo_df <- rbind(parents, children, stringsAsFactors = F)
hpo_df <- as.data.table(hpo_df[!duplicated(hpo_df ),])


index <- data.frame(HPO_ID = unique(hpo_names$HPO_ID), relative = unique(hpo_names$HPO_ID), stringsAsFactors = F)
index$relation <- "index"

#combine into one df
hpo_df <- rbind(hpo_df, index)
hpo_df <-hpo_df[!duplicated(hpo_df ),]
#hpo_df <- merge(hpo_df, hpo_names, by = "HPO_ID")
hpo_df  <- hpo_df [HPO_ID %in% unique(pat_hpo$HPO_ID)]
#hpo_df <- hpo_df[relation != "parent"]


#load gene - hpo 
hpo_gene <- fread("http://compbio.charite.de/jenkins/job/hpo.annotations/lastSuccessfulBuild/artifact/util/annotation/phenotype_to_genes.txt")
colnames(hpo_gene) <- c("HPO_ID", "HPO_name", "entrezID", "geneID", "Additional_Info",  "source", "disease-ID")
hpo_gene <- hpo_gene[, c("HPO_ID","geneID")]
hpo_gene <- hpo_gene[!duplicated(hpo_gene)]


# Annotate relatives by gene
hpo_df <- merge(hpo_df, hpo_gene, by.x= "relative", by.y = "HPO_ID")


# Consider genes matched to relatives - matched to index
hpo_df$relative <- NULL 
hpo_df$relation <- NULL 
hpo_df <- hpo_df[!duplicated(hpo_df), ]
 

pat_hpo <- merge(pat_hpo, hpo_df,  by = "HPO_ID", allow.cartesian=TRUE )
pat_hpo$HPO_ID <- NULL
pat_hpo <-pat_hpo[!duplicated(pat_hpo ),]
pat_hpo[ , HPO_match := T]


#' number of genes matching phenotype per sample
os <- pat_hpo[, .N, by=  SAMPLE_ID]
os <- os[order(N),]
os$rank <- seq(1: nrow(os))
ggplot(os)+
  geom_line( aes( x= rank, y = N), size=1.7)+
  theme_classic()+
  scale_y_continuous(trans='log2')+
  scale_x_continuous(breaks= c(1, 25,  50,  75,  100,  130),  limits=c(1, nrow(os)))+ #
  geom_hline(yintercept = median(os$N), linetype = "dashed")+
  xlab("Sample rank") + 
  ylab("#Genes matching HPO phenotype")




write_tsv(pat_hpo,  snakemake@output$hpo_match)

