#'---
#' title: Supplementary Figure 1e gene set coverage by GTeX proteomics
#' author: smirnovd
#' wb:
#'  input: 
#'  - config: 'src/config.R'
#'  - gencode_annotation: '`sm config["DATASETS"] + "/gene_annotation_v29.tsv"`'
#'  - disease_genes: '`sm config["DATASETS"] + "/disease_genes.tsv"`'
#'  - datasets: '`sm config["DATASETS"]`'
#'  output:
#'  - fig: '`sm config["FIGURE_DIR"] + "/Supplementary_figures/S_Fig1_e.pdf"`'
#'  - detected_proteins_gtex: '`sm config["PROC_DATA"] + "/integration/detected_proteins_gtex.tsv"`'
#' output: 
#'   html_document:
#'    code_folding: hide
#'    code_download: TRUE
#'---

# load config
source(snakemake@input$config)


# Load disease genes table
# dis_genes <- fread('/s/project/mitoMultiOmics/multiOMICs_integration/datasets/disease_genes.tsv')
dis_genes <- fread(snakemake@input$disease_genes)



# Get all protein coding genes
# genecode_v29 <- fread('/s/project/mitoMultiOmics/multiOMICs_integration/datasets/gene_annotation_v29.tsv')
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

rm(protein_coding)



####################################################################


# Download data
# download.file(destfile= paste0("/s/project/mitoMultiOmics/multiOMICs_integration/datasets", "/Table_S1_gene_info_at_protein_level.xlsx" ), url= 'https://storage.googleapis.com/gtex_egtex/proteomics/Table_S1_gene_info_at_protein_level.xlsx')
download.file(destfile= paste0(snakemake@input$datasets, "/Table_S1_gene_info_at_protein_level.xlsx" ), url= 'https://storage.googleapis.com/gtex_egtex/proteomics/Table_S1_gene_info_at_protein_level.xlsx')


# Load GTEx annotation
# sa_gtex <- read_excel(paste0("/s/project/mitoMultiOmics/multiOMICs_integration/datasets", "/Table_S1_gene_info_at_protein_level.xlsx" ), sheet = 2) %>% as.data.table()
sa_gtex <- read_excel(paste0(snakemake@input$datasets, "/Table_S1_gene_info_at_protein_level.xlsx" ), sheet = 2) %>% as.data.table()
sa_gtex[, proteome_ID:= paste0("Run", Run, "_Tag", Tag)]


# Load GTEx proteomics data
# gtex <- read_excel(paste0("/s/project/mitoMultiOmics/multiOMICs_integration/datasets", "/Table_S1_gene_info_at_protein_level.xlsx" ), sheet = 6)
gtex <- read_excel(paste0(snakemake@input$datasets, "/Table_S1_gene_info_at_protein_level.xlsx" ), sheet = 6)
gtex <- gtex[-c(1, 2), -1] 
colnames(gtex)[1] <- "gene_id"
gtex <- merge(genecode_v29[, c("geneID", "gene_id")], gtex, by= "gene_id" )
gtex$gene_id <- NULL

# Convert to numeric
gtex <- cbind(gtex[, 1], sapply(gtex[ , 2 : ncol(gtex)], function(x) {as.numeric(x) }) )
gtex[is.na(gtex)] <- 0
gtex <- rowsum(gtex[,2:ncol(gtex)], gtex$geneID)
gtex <- gtex[!duplicated(gtex), ]
sa_gtex <- sa_gtex[ proteome_ID %in% colnames(gtex)]


# Aggregate some tissues 

sa_gtex[Sample %in% c("Artery - Aorta", "Artery - Coronary", "Artery - Tibial"), Sample := "Artery"]
sa_gtex[Sample %in% c("Brain - Cerebellum" , "Brain - Cortex"), Sample := "Brain"]
sa_gtex[Sample %in% c("Breast - Mammary Tissue" ), Sample := "Breast"]
sa_gtex[Sample %in% c("Colon - Sigmoid"  , "Colon - Transverse" ), Sample := "Colon"]
sa_gtex[Sample %in% c("Esophagus - Gastroesophageal Junction" , "Esophagus - Mucosa", "Esophagus - Muscularis"   ), Sample := "Esophagus"]
sa_gtex[Sample %in% c("Heart - Atrial Appendage" , "Heart - Left Ventricle"  ), Sample := "Heart"]
sa_gtex[Sample %in% c("Muscle - Skeletal"  ), Sample := "Muscle"]
sa_gtex[Sample %in% c("Nerve - Tibial"  ), Sample := "Nerve"]
sa_gtex[Sample %in% c("Skin - Not Sun Exposed (Suprapubic)", "Skin - Sun Exposed (Lower leg)" ), Sample := "Skin"]
sa_gtex[Sample %in% c("Small Intestine - Terminal Ileum" ), Sample := "Small Intestine"]
sa_gtex[Sample %in% c("Minor Salivary Gland" ), Sample := "Salivary Gland"]

##########################################
gtex_all <- data.frame()
detected_gtex <- data.frame()

for (tissue in unique(sa_gtex[Sample != "reference"]$Sample)){
  tissue_samples <- sa_gtex[Sample ==  tissue ]$proteome_ID
  expr <- gtex[ , tissue_samples]
  
  expr$geneID <-  rownames(expr)
  expr_p = melt(expr, id.vars=c("geneID")) %>% as.data.table()
  expr_p <- expr_p[value > 0 , .N, by = "geneID"]
  expr_p[, half :=  N >= max(N)/2]
  expr_protein <- expr_p
  expr_p$N <- NULL
  
  expr_protein[ , once :=  N == 1]
  expr_protein[ , all :=  N == max(N)]
  expr_protein[ , TISSUE :=  tissue]
  detected_gtex <-rbind(detected_gtex , expr_protein)

  

  dg_expr <- merge(dis_genes, expr_p, by = "geneID")
  dg_expr[, TISSUE := tissue]
  dg_expr[ half == T, HALF := .N, by = DISEASE ]
  dg_expr <- dg_expr[, c("TISSUE", "DISEASE", "total", "HALF" )]
  dg_expr <- dg_expr[!duplicated(dg_expr)]
  dg_expr <- dg_expr[!is.na(HALF)]
  
  detected_expr <- melt(dg_expr, id.vars=c("TISSUE", "DISEASE", "total" )) %>% as.data.table()
  detected_expr <- detected_expr[!duplicated(detected_expr)]
  setnames(detected_expr, c("variable", "value" ), c("DETECTED", "N"))
  detected_expr[, prop:= N / total ]
  gtex_all <-rbind(gtex_all , detected_expr)
}


# Write detected proteins
# write_tsv(detected_gtex, '/s/project/mitoMultiOmics/multiOMICs_integration/processed_data/integration/detected_proteins_gtex.tsv')
write_tsv(detected_gtex,  snakemake@output$detected_proteins_gtex)




gtex_all$dis_n <- paste0(gtex_all$DISEASE, "\n"  ,"(",gtex_all$total, ")")
# unique(gtex_all$dis_n)
gtex_all$dis_n <- factor(gtex_all$dis_n, levels = c("Protein coding\n(20336)", 
                                                    "MITO\n(413)", 
                                                    "Neuromuscular\n(132)", 
                                                    "Neurology\n(284)", 
                                                    "Hematology\n(50)",
                                                    "Ophthalmology\n(366)",
                                                    "OMIM\n(4354)"))

gtex_all[ , CAT:= "No" ]
gtex_all[TISSUE %in% c("Muscle", "Skin") , CAT:= "Yes" ]






#+ fig.width=10, fig.height=6
s_fig <- ggplot(gtex_all, aes(dis_n, prop)) + 
  geom_boxplot( ) + # fill = "#FB9A99" 
  geom_point( aes(color = TISSUE, shape = CAT) , position = "jitter", size = 2) + 
  scale_y_continuous( labels=scales::percent) +
  scale_fill_brewer() +
  labs( y = "Median coverage")+
  ggtitle("Protein detection across all GTEx tissues")+
  theme_classic()+
  scale_shape_manual( values = c(19, 15))+
  theme(legend.position="bottom", 
        plot.title = element_text(hjust = 0.5, size=12,face="bold"),
        axis.title.y = element_text(face="bold", size=12) , 
        axis.title.x = element_blank() ,
        axis.text.x = element_text(size=12, face="bold") ,
        legend.direction = "horizontal", 
        axis.text.y = element_text(face="bold", size=12, hjust = 0.5))


#s_fig

pdf(snakemake@output$fig, # "/s/project/mitoMultiOmics/multiOMICs_integration/Figures/Supplementary_figures/S_Fig1_e.pdf",  
    width = 10, height = 6,  useDingbats=FALSE )
print(s_fig) 
dev.off()

