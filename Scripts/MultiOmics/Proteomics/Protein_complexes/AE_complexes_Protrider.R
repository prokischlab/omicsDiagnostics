#'---
#' title: Aberrantly expressed CORUM complexes PROTRIDER
#' author: Dmitrii Smirnov
#' wb:
#'  input:
#'  - config: 'src/config.R'
#'  - protrider: '`sm config["PROC_DATA"] + "/protrider/PROTRIDER_results.rds"`'
#'  - sample_annotation: '`sm config["ANNOTATION"]`'
#'  - corum: '`sm config["PROC_DATA"] + "/Complexes/CORUM.tsv"`'
#'  output:
#'  - complex_results: '`sm config["PROC_DATA"] + "/Complexes/Complex_outliers_PROTRIDER.rds"`'
#'  - subunits_result: '`sm config["PROC_DATA"] + "/Complexes/aberrant_complex_subunits_PROTRIDER.tsv"`'
#' output:
#'   html_document:
#'    code_folding: hide
#'    code_download: TRUE
#'---


source(snakemake@input$config)
#load functions
source("src/functions/LIMMA/limma_functions.R")

# READ ANNOTATION
sa <- fread(snakemake@input$sample_annotation)
sa <- sa[USE_FOR_PROTEOMICS_PAPER == T]

# Read protein outlier results
prot <- readRDS(snakemake@input$protrider) %>% as.data.table()


# Subset necessary columns 
prot <- prot[,  c("SAMPLE_ID", "geneID","PROTEIN_LOG2FC")]

# Read CORUM complexes
corum <- fread(snakemake@input$corum)

paste("Number of unique protein complexes:", uniqueN(corum$COMPLEX) )
# Keep only proteins with known PC
paste( "proteins analysed by AE", uniqueN(prot$geneID))
prot <- prot[geneID %in% unique(corum$geneID)]
paste( "proteins analysed by AE, detected in complexes", uniqueN(prot$geneID))

# Merge and remove duplicates
prot_pc <- merge(corum, prot, by = "geneID", allow.cartesian=TRUE)
prot_pc <- prot_pc[ !duplicated(prot_pc) ]

# Count subunits detected by proteomics
prot_pc[, N_quantified_subunits:= .N, by = .(COMPLEX, SAMPLE_ID)]

# Use complexes with at least two quantified subunits
prot_pc <- prot_pc[ N_quantified_subunits >=2 ]

# Use complexes with at least half of subunits detected
prot_pc <- prot_pc[ N_quantified_subunits >= N_subunits/2 ]


paste("Number of unique protein complexes after filterring:", uniqueN(prot_pc$COMPLEX) )
paste( "Number of proteins in PC analysis", uniqueN(prot_pc$geneID))


# Compute mean FC per protein complex, per sample
prot_pc[, mean_COMPLEX_LOG2FC:= mean(PROTEIN_LOG2FC, na.rm = T), by = .(SAMPLE_ID, COMPLEX)]


complex <- prot_pc[ , c("SAMPLE_ID", "COMPLEX", "N_subunits", "N_quantified_subunits", "mean_COMPLEX_LOG2FC" )]
complex <- complex[ !duplicated(complex) ]

# Mean fold change is normally distributed
hist(complex$mean_COMPLEX_LOG2FC, breaks = 80)

# Example distribution of mean fold change per complex
hist(complex[COMPLEX == "Respiratory chain complex I (holoenzyme), mitochondrial"]$mean_COMPLEX_LOG2FC)


# Compute mean, standard deviation, FC, Z-score and p-values per complex by fitting normal distribution
complex[ , MEAN := fitdistr(na.exclude(mean_COMPLEX_LOG2FC), "normal")$estimate[1] , by = COMPLEX]
complex[ , SD := fitdistr(na.exclude(mean_COMPLEX_LOG2FC), "normal")$estimate[2] , by = COMPLEX ]
complex[ , COMPLEX_LOG2FC := mean_COMPLEX_LOG2FC - MEAN ] 
complex[ , COMPLEX_ZSCORE := COMPLEX_LOG2FC/SD ]
complex[ , COMPLEX_FC := 2^COMPLEX_LOG2FC ] 
complex[ , Pval := pnorm(mean_COMPLEX_LOG2FC , mean = MEAN, sd = SD ) ] 

# convert to two-tailed p-values
complex[ , COMPLEX_PVALUE := 2*pmin(Pval, 1-Pval) ] 
hist(complex$COMPLEX_PVALUE)
complex$Pval <- NULL

# Define outliers by p-value or z-score
complex[ , COMPLEX_PADJ := p.adjust(COMPLEX_PVALUE, method = 'BY'), by = SAMPLE_ID ]
complex[, COMPLEX_outlier :=  COMPLEX_PADJ < 0.1 ]  
complex[, COMPLEX_Z_outlier :=   abs(COMPLEX_ZSCORE) >=3 ]  
complex$mean_COMPLEX_LOG2FC <- NULL



#' ### Number of samples with at least one aberrant complex 
paste("Total samples:", complex[, uniqueN(SAMPLE_ID)])
paste("Significance based:", complex[COMPLEX_outlier == T, uniqueN(SAMPLE_ID)])
paste("Z-score based:", complex[COMPLEX_Z_outlier == T, uniqueN(SAMPLE_ID)])


#' ## Aberrant complexes per sample
plotAberrantProteinPerSample(complex[ COMPLEX_outlier == T , .N, by = SAMPLE_ID])
plotAberrantProteinPerSample(complex[ COMPLEX_Z_outlier == T , .N, by = SAMPLE_ID])


# Annotate with complex description 
corum_annotate <- fread("curl http://mips.helmholtz-muenchen.de/corum/download/allComplexes.txt.zip | funzip")
corum_annotate <- corum_annotate[Organism == 'Human']
corum_annotate <- corum_annotate[ComplexName %in% unique(complex$COMPLEX), 
               c("ComplexName", "Complex comment", "Disease comment",
                   "GO description", "FunCat description",
                   "Synonyms" ,"subunits(Gene name)" )]
corum_annotate <- corum_annotate[ !duplicated(corum_annotate) ]
corum_annotate <- corum_annotate[ !duplicated(corum_annotate$ComplexName) ]

res <- merge(complex, corum_annotate, by.x= "COMPLEX", by.y = "ComplexName", all.x = T)
res <- res[order(SAMPLE_ID)]

saveRDS(res, snakemake@output$complex_results)


#' ## Results for aberrantly expressed CORUM complexes
DT::datatable(res[ COMPLEX_outlier == T | COMPLEX_Z_outlier == T], 
              caption = "Aberrantly expressed CORUM complexes", style = 'bootstrap', filter = 'top', escape = F,
              extensions = c( 'Buttons', 'ColReorder' ), 
              options = list( colReorder = TRUE, dom = 'Bfrtip',
                              buttons = c('copy', 'csv', 'excel', 'pdf', 'print')))  




# Annotate proteins with complex outlier evidence
corum$N_subunits <- NULL
subunits <- merge(corum, complex, by = "COMPLEX", all =T,  allow.cartesian=TRUE)
subunits <- subunits[!duplicated(subunits)]

#' ## Aberrant PC subunits per sample
plotAberrantProteinPerSample(subunits[ COMPLEX_outlier == T , .N, by = c('SAMPLE_ID')])
plotAberrantProteinPerSample(subunits[ COMPLEX_Z_outlier == T , .N, by = c('SAMPLE_ID')])


# Subset only significant
subunits <- subunits[ COMPLEX_outlier == T | COMPLEX_Z_outlier == T]


# Annotate with causal genes 
sa$sample_gene <- paste0(sa$SAMPLE_ID, "_", sa$KNOWN_MUTATION)
subunits$sample_gene <- paste0(subunits$SAMPLE_ID, "_", subunits$geneID)
subunits[, causal_gene := sample_gene %in% unique(sa$sample_gene)]


# Results for subunits aberrantly expressed CORUM complexes
DT::datatable(subunits,
              caption = "Aberrantly expressed CORUM complex subunits", style = 'bootstrap', filter = 'top', escape = F,
              extensions = c( 'Buttons', 'ColReorder' ),
              options = list( colReorder = TRUE, dom = 'Bfrtip',
                              buttons = c('copy', 'csv', 'excel', 'pdf', 'print')))


subunits <- subunits[, .(sample_gene, COMPLEX)]
subunits <- subunits[!duplicated(subunits)]
subunit <- aggregate(subunits[, -1], by= list(subunits$sample_gene), paste)
setnames(subunit, c("sample_gene" , "Aberrant_complexes") )
subunit$Aberrant_complexes <- as.character(subunit$Aberrant_complexes)


write_tsv(subunit,  snakemake@output$subunits_result)

