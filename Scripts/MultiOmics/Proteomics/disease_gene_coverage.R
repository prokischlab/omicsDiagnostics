#'---
#' title: Disease gene coverage by proteomics
#' author: Dmitrii Smirnov
#' wb:
#'  input: 
#'  - gencode_annotation: '`sm config["DATASETS"] + "/gene_annotation_v29.tsv"`'
#'  - disease_genes: '`sm config["DATASETS"] + "/disease_genes.tsv"`'
#'  - sample_annotation: '`sm config["ANNOTATION"]`'
#'  - norm_data: '`sm config["PROC_DATA"] + "/limma/proteomics_normalized_not_imputed.tsv"`'
#'  output:
#'  - detected_proteins: '`sm config["PROC_DATA"] + "/integration/detected_proteins.tsv"`'
#'  - protein_coverage: '`sm config["PROC_DATA"] + "/integration/protein_coverage.tsv"`'
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



# Read protein matrix 
proteomics <- fread(snakemake@input$norm_data)
colnames(proteomics)[1] <-  "geneID"

proteomics_detected = melt(proteomics, id.vars= "geneID") %>% as.data.table()

# Subset for paper cases
proteomics_detected <- proteomics_detected[variable %in% sa$SAMPLE_ID ]

# Count genes detected once ...
proteomics_detected <- proteomics_detected[value > 0 , .N, by = "geneID"]
proteomics_detected[, once := T]
proteomics_detected[, half :=  N >= max(N)/2]
proteomics_detected[, all :=  N == max(N)]
proteomics_detected$N <- NULL


write_tsv(proteomics_detected,  snakemake@output$detected_proteins)


# Combine with dis gene list
dis_genes_fib <- merge(dis_genes, proteomics_detected, by = "geneID")
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


#' # Proteomics coverage in fibroblasts
#+ fig.width=7, fig.height=4
ggplot(detected_fib, aes(dis_n, prop)) + 
  geom_col(aes(fill = DETECTED)) +  # stat= 'identity', 
  scale_y_continuous( labels=scales::percent) +
  scale_fill_brewer(palette="Reds") +
  labs( y = "fraction covered")+
  theme_classic()+
  theme(legend.position="none",  
        axis.title.y = element_text(face="bold", size=12) , 
        axis.title.x = element_blank() ,
        axis.text.x = element_text(size=12, face="bold") ,
        legend.title = element_blank(), legend.direction = "horizontal", 
        axis.text.y = element_text(face="bold", size=12, hjust = 0.5),
        plot.margin = margin(0, 0, 0, 0, "cm")) 


write_tsv(detected_fib,  snakemake@output$protein_coverage)

