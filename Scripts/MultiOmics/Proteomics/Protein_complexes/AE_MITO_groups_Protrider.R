#'---
#' title: Aberrantly expressed HGNC mito groups PROTRIDER
#' author: Dmitrii Smirnov
#' wb:
#'  input:
#'  - config: 'src/config.R'
#'  - protrider: '`sm config["PROC_DATA"] + "/protrider/PROTRIDER_results.rds"`'
#'  - sample_annotation: '`sm config["ANNOTATION"]`'
#'  - mito_groups: '`sm config["DATASETS"] + "/HGNC_mito_groups.tsv"`'
#'  output:
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

# Read HGNC mito groups
compl <- fread(snakemake@input$mito_groups)

setnames(compl, c("Group_name", "COMPLEX"), c("COMPLEX", "COMPLEX2") )
compl[, N_subunits := .N, by = COMPLEX]
compl <- compl[N_subunits >1 ]

paste("Number of unique Mito complexes/groups:", uniqueN(compl$COMPLEX) )
# Keep only proteins with known PC
paste( "proteins analysed by AE", uniqueN(prot$geneID))
prot <- prot[geneID %in% unique(compl$geneID)]
paste( "proteins analysed by AE, detected in Mito complexes/groups", uniqueN(prot$geneID))

# Merge and remove duplicates
prot_pc <- merge(compl, prot, by = "geneID", allow.cartesian=TRUE)
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


complex <- prot_pc[ , c("SAMPLE_ID", "COMPLEX","COMPLEX2", "N_subunits", "N_quantified_subunits", "mean_COMPLEX_LOG2FC" )]
complex <- complex[ !duplicated(complex) ]

# Mean fold change is normally distributed
hist(complex$mean_COMPLEX_LOG2FC, breaks = 80)

# Example distribution of mean fold change per complex
hist(complex[COMPLEX == "Mitochondrial complex III: ubiquinol-cytochrome c reductase complex subunits" ]$mean_COMPLEX_LOG2FC)
unique(complex$COMPLEX)

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
#'  P-value < 0.1
plotAberrantProteinPerSample(complex[ COMPLEX_outlier == T , .N, by = SAMPLE_ID])
#' |Z-score| >= 3
plotAberrantProteinPerSample(complex[ COMPLEX_Z_outlier == T , .N, by = SAMPLE_ID])


complex <- complex[order(SAMPLE_ID)]


#' # Results for aberrantly expressed HGNC mito groups
DT::datatable(complex[ COMPLEX_outlier == T | COMPLEX_Z_outlier == T], 
              caption = "Aberrantly expressed HGNC mito groups", style = 'bootstrap', filter = 'top', escape = F,
              extensions = c( 'Buttons', 'ColReorder' ), 
              options = list( colReorder = TRUE, dom = 'Bfrtip',
                              buttons = c('copy', 'csv', 'excel', 'pdf', 'print')))  



#' # Plots per case 
#' Nominal p-values inicated in the figures
sax <- sa[!is.na(KNOWN_MUTATION) & !is.na(CATEGORY)]
sax <- sax[order(CATEGORY)]
complexX <- complex[COMPLEX2 %in% c("mRibo", "RCCI_assembly", "RCCI", "RCCII", "RCCIII", "RCCIV", "RCCV")]
sax <- sax[SAMPLE_ID %in% unique(complexX$SAMPLE_ID)]

#+ fig.width=8, fig.height=5
for (i in unique(sax$SAMPLE_ID)){
  complexX[ , case := F]
  complexX[ SAMPLE_ID == i, case := T]
  complexX[ , Pval := round(COMPLEX_PVALUE,  3)]
  print(paste("Sample ID:", i, "Causal gene:", sax[SAMPLE_ID == i]$KNOWN_MUTATION ) )
  
  plot_compl <- ggplot(complexX, aes(COMPLEX2,  COMPLEX_ZSCORE ))+
    geom_hline(yintercept = 0, color = "grey50") +
    geom_quasirandom(data = complexX[case == F ] , aes(COMPLEX2,  COMPLEX_ZSCORE),color = "grey80" , size=0.5) +
    stat_summary(fun.data = mean_sdl, fun.args = list(mult=1), geom="errorbar", color="black", width=0.2) +
    stat_summary(fun = mean, geom="point", color="black") +
    ylab("Protein complex zScore") + 
    geom_point( data = complexX[case == T ], aes(COMPLEX2,  COMPLEX_ZSCORE), colour = "darkorange" , size = 2) +
    geom_text_repel( data = complexX[case == T ], aes(COMPLEX2,  COMPLEX_ZSCORE, label = Pval) ) +
    theme_classic()+
    theme( legend.title = element_blank(),
           axis.title.x = element_blank(),
           plot.title = element_text(face="bold",  size=14, hjust = 0.5),
           axis.title.y = element_text(face="bold",  size=10),
           axis.text.x = element_text(face="bold",  size=9, angle = 45, hjust = 1),
           axis.text.y = element_text(face="bold",  size=9),
           legend.position = "none")
  print(plot_compl)
  
}







