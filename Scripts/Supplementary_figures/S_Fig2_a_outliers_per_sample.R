#'---
#' title: Supplementary Fig 2a outliers per sample 
#' author: smirnovd
#' wb:
#'  input: 
#'  - config: 'src/config.R'
#'  - sample_annotation: '`sm config["ANNOTATION"]`'
#'  - outrider_results: '`sm config["PROC_DATA"] + "/outrider/OUTRIDER_results.rds"`'
#'  - protrider_results: '`sm config["PROC_DATA"] + "/protrider/PROTRIDER_results.rds"`'
#'  - limma_results: '`sm config["PROC_DATA"] + "/limma/LIMMA_results.rds"`'
#'  output:
#'  - fig: '`sm config["FIGURE_DIR"] + "/Supplementary_figures/S_Fig2_a.pdf"`'
#' output: 
#'   html_document:
#'    code_folding: hide
#'    code_download: TRUE
#'---

# load config
source(snakemake@input$config)


# sa <- fread('/s/project/mitoMultiOmics/multiOMICs_integration/raw_data/proteomics_annotation.tsv')
# limma <- readRDS('/s/project/mitoMultiOmics/multiOMICs_integration/processed_data/limma/LIMMA_results.rds') %>% as.data.table()
# protr <- readRDS("/s/project/mitoMultiOmics/multiOMICs_integration/processed_data/protrider/PROTRIDER_results.rds") %>% as.data.table()
# outrider <- readRDS("/s/project/mitoMultiOmics/multiOMICs_integration/processed_data/outrider/OUTRIDER_results.rds") %>% as.data.table()

# Load sample annotation
sa <- fread(snakemake@input$sample_annotation)
sa <- sa[USE_FOR_PROTEOMICS_PAPER == T]
# Keep only full sample list
sa <- data.frame(sa[, "SAMPLE_ID"])

# Create outliers per sample df
os <- data.frame() 

##################
# Load LIMMA results
##################

limma <- readRDS(snakemake@input$limma_results)  %>% as.data.table()
limma <- limma[!is.na(PROTEIN_PVALUE)]
paste("Number of outliers in LIMMA:" ,nrow(limma[ PROTEIN_outlier == T]))

# Keep only cases with outliers and count by sample
osx <- limma[ PROTEIN_outlier == T , .N, by = c('SAMPLE_ID')]
osx <- as.data.table(merge(sa, osx, by = "SAMPLE_ID", all.x = T))
# Set N outliers as 0 for cases wo outliers
osx[is.na(N), N:= 0]
osx <- osx[order(N),]
osx$Method <- rep("LIMMA", nrow(osx) )
osx$Rank <- seq( 1: nrow(osx) )
os <- rbind(os, osx)


##################
#Load PROTRIDER results
##################
protr <- readRDS(snakemake@input$protrider_results) %>% as.data.table()

protr <- protr[!is.na(PROTEIN_PVALUE)]
paste("Number of outliers in PROTRIDER:", nrow(protr[ PROTEIN_outlier == T]) )

osx <- protr[ PROTEIN_outlier == T , .N, by = c('SAMPLE_ID')]
osx <- as.data.table(merge(sa, osx, by = "SAMPLE_ID", all.x = T))
osx[is.na(N), N:= 0]
osx <- osx[order(N),]
osx$Method <- rep("PROTRIDER", nrow(osx) )
osx$Rank <- seq( 1: nrow(osx) )
os <- rbind(os, osx)


##################
#Load OUTRIDER results
##################


outrider <- readRDS(snakemake@input$outrider_results) %>% as.data.table()

outrider <- outrider[!is.na(RNA_PVALUE)]
paste("Number of outliers in OUTRIDER:", nrow(outrider[ RNA_outlier == T]) )


osx <- outrider[ RNA_outlier == T , .N, by = c('SAMPLE_ID')]
osx <- as.data.table(merge(sa, osx, by = "SAMPLE_ID", all.x = T))
osx[is.na(N), N:= 0]
osx <- osx[order(N),]
osx$Method <- rep("OUTRIDER", nrow(osx) )
osx$Rank <- seq( 1: nrow(osx) )
os <- rbind(os, osx)



###########################################################
s_fig <- ggplot(os, aes(Rank , N, color = Method))+
  geom_line( size = 1.2)+
  theme_classic()+
  scale_color_manual( values = c("OUTRIDER" = "#A6CEE3",
                                 "PROTRIDER" = "#FB9A99",
                                 "LIMMA" = "darkorange"))+

  geom_hline(yintercept = round(median(os[Method == "LIMMA"]$N)), color = "darkorange", linetype = "dashed")+
  geom_hline(yintercept = round(median(os[Method == "PROTRIDER"]$N)), color = "#FB9A99", linetype = "dashed")+
  geom_hline(yintercept = round(median(os[Method == "OUTRIDER"]$N)), color = "#A6CEE3", linetype = "dashed")+

  scale_y_continuous(breaks= c(2, 
                               round(median(os[Method == "OUTRIDER"]$N)), 
                               round(median(os[Method == "LIMMA"]$N)),
                               round(median(os[Method == "PROTRIDER"]$N)), 
                               10,
                               25,
                               50, 
                               100,
                               max(os$N)),  trans='log2')+
  scale_x_continuous(breaks= c(1, 25,  50,  75,  100,  125,  147), limits=c(1, nrow(sa)))+
  xlab("Sample rank") + 
  ylab("# aberrantly expressed genes")+
  theme(plot.title = element_text(hjust = 0.5),
        axis.title.x= element_text( size=10, margin = NULL,face="bold"),
        axis.title.y= element_text( size=10, margin = NULL,face="bold"),
        axis.text.x = element_text(face="bold",  size=10),
        axis.text.y = element_text(face="bold",  size=10),
        legend.title = element_blank(), # legend.position = c(0.2, 0.8)
        legend.position = "left")

#+ fig.width=6, fig.height=6
s_fig


pdf(snakemake@output$fig, #"/s/project/mitoMultiOmics/multiOMICs_integration/Figures/Supplementary_figures/S_Fig2_a.pdf",  
    width = 8, height =5,  useDingbats=FALSE )
print(s_fig) 
dev.off()





