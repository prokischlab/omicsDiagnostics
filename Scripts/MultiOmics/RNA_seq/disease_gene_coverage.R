#'---
#' title: Disease gene coverage by RNAseq
#' author: Dmitrii Smirnov
#' wb:
#'  input: 
#'  - gencode_annotation: '`sm config["DATASETS"] + "/gene_annotation_v29.tsv"`'
#'  - disease_genes: '`sm config["DATASETS"] + "/disease_genes.tsv"`'
#'  - sample_annotation: '`sm config["ANNOTATION"]`'
#'  - ods_unfiltered: '`sm config["PROC_DATA"] + "/outrider/ods_unfiltered.Rds"`'
#'  output:
#'  - detected_transcripts: '`sm config["PROC_DATA"] + "/integration/detected_transcripts.tsv"`'
#'  - rna_coverage: '`sm config["PROC_DATA"] + "/integration/rna_coverage.tsv"`'
#' output: 
#'   html_document:
#'    code_folding: hide
#'    code_download: TRUE
#'---

# Load config
source('src/config.R')


# Load disease genes table
dis_genes <- fread(snakemake@input$disease_genes)


# Get all protein coding genes
genecode_v29 <- fread(snakemake@input$gencode_annotation)
genecode_v29[, geneID := toupper(gene_name_unique)]
protein_coding <- genecode_v29[gene_type == 'protein_coding']
protein_coding[, DISEASE := 'Protein coding']
protein_coding <- protein_coding[geneID != "" & !is.na(geneID) , c("geneID", "DISEASE")]
protein_coding <- protein_coding[!duplicated(protein_coding)]
protein_coding[, total := .N]
protein_coding[, ORIGIN := 'genecode v29']

# Combine
dis_genes <- rbind( dis_genes, protein_coding)
dis_genes[ , geneID := toupper(geneID)]
dis_genes <- dis_genes[!duplicated(dis_genes)]

# Subset only protein coding genes
dis_genes <- dis_genes[geneID %in% unique(protein_coding$geneID)]
dis_genes[ , total := .N, by = DISEASE]
rm(protein_coding, genecode_v29)

dis_genes <- dis_genes[DISEASE %in% c("Protein coding", "MITO", "Neuromuscular", "Neurology", "OMIM" ) ]

####################################################################


# READ ANNOTATION
sa <- fread(snakemake@input$sample_annotation)
sa <- sa[USE_FOR_PROTEOMICS_PAPER == T]


### RNA unfiltered
ods <- readRDS(snakemake@input$ods_unfiltered)
ods_100 <- OUTRIDER::filterExpression(ods, percentile = 0)   # detected in all
ods_50 <- OUTRIDER::filterExpression(ods, percentile = 0.5)  # detected in 50%
ods_once <- OUTRIDER::filterExpression(ods, percentile = 0.9999) # detected in 0,001% (once)
rm(ods)

rnaseq <- as.data.table(rowData(ods_once))
rnaseq[, once := T]
rnaseq <- rnaseq[ , c("gene_name_unique", "once")]
rnaseq[ , half := gene_name_unique %in% rowData(ods_50)$gene_name_unique ]
rnaseq[ , all := gene_name_unique %in% rowData(ods_100)$gene_name_unique ]
setnames(rnaseq, "gene_name_unique", "geneID")
rnaseq[ , geneID := toupper(geneID)]
rnaseq <- rnaseq[!duplicated(rnaseq)]


write_tsv(rnaseq,  snakemake@output$detected_transcripts)




# Combine with dis gene list
rnaseq <- rnaseq[geneID %in% unique(dis_genes$geneID) ]
dis_genes_fib <- merge(dis_genes, rnaseq, by = "geneID")
dis_genes_fib[, TISSUE := "Fibroblasts"]
dis_genes_fib[ once == T, ONCE_1 := .N, by = DISEASE ]
dis_genes_fib[ half == T, HALF_1 := .N, by = DISEASE ]
dis_genes_fib[ all == T, ALL := .N, by = DISEASE ]
dis_genes_fib <- dis_genes_fib[, c("TISSUE", "DISEASE", "total", "ONCE_1","HALF_1","ALL"  )]
dis_genes_fib <- dis_genes_fib[!duplicated(dis_genes_fib)]
dis_genes_fib[, HALF := HALF_1 - ALL]
dis_genes_fib[, ONCE := ONCE_1 - HALF_1]
dis_genes_fib <- dis_genes_fib[, c("TISSUE", "DISEASE", "total", "ONCE","HALF","ALL"  )]
dis_genes_fib <- dis_genes_fib[!is.na(ALL)]

detected_fib <- melt(dis_genes_fib, id.vars=c("TISSUE", "DISEASE", "total" )) %>% as.data.table()
detected_fib <- detected_fib[!duplicated(detected_fib)]
setnames(detected_fib, c("variable", "value" ), c("DETECTED", "N"))
detected_fib[, prop:= N / total ]


detected_fib$DISEASE <- factor(detected_fib$DISEASE, levels = c("Protein coding", "MITO", "Neuromuscular", "Neurology", "OMIM" ))

detected_fib$dis_n <- paste0(detected_fib$DISEASE, "\n"  ,"(",detected_fib$total, ")")
# unique(detected_fib$dis_n)
detected_fib$dis_n <- factor(detected_fib$dis_n, levels = c("Protein coding\n(20336)", "MITO\n(388)", "Neuromuscular\n(132)", "Neurology\n(284)", "OMIM\n(4270)" ))



#' # RNA-seq coverage in fibroblasts
#+ fig.width=7, fig.height=4
ggplot(detected_fib, aes(dis_n, prop)) + 
  geom_col(aes(fill = DETECTED)) +  # stat= 'identity', 
  scale_y_continuous( labels=scales::percent) +
  scale_fill_brewer(palette="Blues") +
  labs( y = "fraction covered")+
  theme_classic()+
  theme(legend.position="none",  
        axis.title.y = element_text(face="bold", size=12) , 
        axis.title.x = element_blank() ,
        axis.text.x = element_text(size=12, face="bold") ,
        legend.title = element_blank(), legend.direction = "horizontal", 
        axis.text.y = element_text(face="bold", size=12, hjust = 0.5),
        plot.margin = margin(0, 0, 0, 0, "cm")) 

write_tsv(detected_fib,  snakemake@output$rna_coverage)


