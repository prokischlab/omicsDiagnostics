#'---
#' title: Supplementary Figure 1d, detected proteins
#' author: Robert Kopajtich, Dmitrii Smirnov
#' wb:
#'  input: 
#'  - config: 'src/config.R'
#'  - raw_prot: '`sm config["RAW_Protein"]`'
#'  - sample_annotation: '`sm config["ANNOTATION"]`'
#'  output:
#'  - fig: '`sm config["FIGURE_DIR"] + "/Supplementary_figures/S_Fig1_d.pdf"`'
#' output: 
#'   html_document:
#'    code_folding: hide
#'    code_download: TRUE
#'---

source(snakemake@input$config)

# Load sample annotation
# sa <- fread('/s/project/mitoMultiOmics/multiOMICs_integration/raw_data/proteomics_annotation.tsv')
sa <- fread(snakemake@input$sample_annotation)
sa <- sa[USE_FOR_PROTEOMICS_PAPER == T]


# Read raw protein matrix
raw_prot <- fread(snakemake@input$raw_prot)
# raw_prot <- fread('/s/project/mitoMultiOmics/multiOMICs_integration/raw_data/proteomics_not_normalized.tsv')
dim(raw_prot)

# Subset for the paper cases
keep_columns <- c( "geneID", sa$SAMPLE_ID)
raw_prot <- raw_prot[ , ..keep_columns ]

prot_det <- as.matrix(raw_prot[,2:ncol(raw_prot), with=FALSE] > 0)




# clean data and sort by num of detected samples
prot_det <- prot_det[,colSums(prot_det) > 5000]     # filter for samples with less than 5000 x TRUE
# prot_det <- prot_det[,order(colSums(prot_det))]   # sort by number of proteins detected


# iterate through the index of you samples 1 ... 200
# and save your last value (here total detected proteins till index i)
all_detected        <- prot_det[,1]
shared_detected     <- prot_det[,1]
cum_num_detected    <- numeric(ncol(prot_det))
shared_num_detected <- numeric(ncol(prot_det))

for(i in 1:ncol(prot_det)){
    all_detected    <- all_detected    | prot_det[,i]
    shared_detected <- shared_detected & prot_det[,i]
    cum_num_detected[i] <- sum(all_detected)
    shared_num_detected[i] <- sum(shared_detected)
}

# merge data and prepare it for ggplot
dt2plot <- melt(1:2, value.name="NumProteins", variable.name="Variable",
                data=data.table(
                    rank=1:ncol(prot_det),
                    sampleID=colnames(prot_det),
                    NumDetect=colSums(prot_det),
                    NumShared=shared_num_detected,
                    NumCum=cum_num_detected
                ))

# plot data
s_fig <- ggplot(dt2plot, aes(rank, NumProteins, color=Variable)) +
    geom_line() +
    ylim(c(0, max(dt2plot[,NumProteins]))) +
    ylab("# of proteins") +
    xlab("# of samples")+
    ggtitle("detected proteins, cummulative, per sample, shared")+
    theme_classic()+
    theme(plot.title = element_text(hjust = 0.5, size=12,face="bold"),
          axis.title.y = element_text(face="bold", size=12) , 
          axis.title.x = element_text(face="bold", size=12) , 
          axis.text.x = element_text(size=12, face="bold") ,
          axis.text.y = element_text(face="bold", size=12),
          legend.title = element_blank(), 
          plot.margin = margin(0, 0, 0, 0, "cm"))

#+ fig.width=6, fig.height=3
s_fig


pdf(snakemake@output$fig, # "/s/project/mitoMultiOmics/multiOMICs_integration/Figures/Supplementary_figures/S_Fig1_d.pdf",  
    width = 3, height =6,  useDingbats=FALSE )
print(s_fig) 
dev.off()



