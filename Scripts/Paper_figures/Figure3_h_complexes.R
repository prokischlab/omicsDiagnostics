#'---
#' title: Figure 3h RCC analysis
#' author: Dmitrii Smirnov
#' wb:
#'  input:
#'  - config: 'src/config.R'
#'  - patient_omics: '`sm config["PROC_DATA"] + "/integration/patient_omics.RDS"`'
#'  - mito_groups: '`sm config["DATASETS"] + "/HGNC_mito_groups.tsv"`'
#'  output:
#'  - fig: '`sm config["FIGURE_DIR"] + "/Fig3_h.pdf"`'
#' output: 
#'   html_document:
#'    code_folding: hide
#'    code_download: TRUE
#'---

source(snakemake@input$config)


# Read integrated omics file 
rp <- readRDS(snakemake@input$patient_omics) %>% as.data.table()


# Read mito complexes 
compl <- fread(snakemake@input$mito_groups)
compl <- compl[ COMPLEX %in% c("mRibo", "RCCI", "RCCII", "RCCIII", "RCCIV", "RCCV"), ]
compl[, gt := "nDNA-encoded"]
compl[Chromosome == "mitochondria", gt := "mtDNA-encoded"]
compl <- compl[ , .( geneID, gt)]

rp <- rp[ geneID %in% compl$geneID, .(SAMPLE_ID, geneID, RNA_FC, PROTEIN_FC)]
rp <- merge(rp, compl, by = "geneID")
rp[ , case := "controls"]
rp[SAMPLE_ID == "OM91786" , case := "subject"]

rp[, mean_RNA_FC:= mean(RNA_FC, na.rm = T), by = .(gt, case)] #SAMPLE_ID, 
rp[, RNA_SD:= sd(PROTEIN_FC, na.rm = T), by = .(gt, case)] #SAMPLE_ID, 

rp[, mean_PROTEIN_FC:= mean(PROTEIN_FC, na.rm = T), by = .(gt, case)] #SAMPLE_ID, 
rp[, PROTEIN_SD:= sd(PROTEIN_FC, na.rm = T), by = .( gt, case)] #SAMPLE_ID, 

rpx <- rp[ , .(gt, case, mean_RNA_FC,RNA_SD,  mean_PROTEIN_FC, PROTEIN_SD)]
rpx <- rpx[!duplicated(rpx)]

rp_r <- rpx[ , .( gt, case, mean_RNA_FC, RNA_SD)]
colnames(rp_r) <- c( "gt", "case", "mean_fold_change", "SD" )
rp_r[ , assay := "RNA"]

rp_p <- rpx[ , .(gt, case, mean_PROTEIN_FC, PROTEIN_SD)]
colnames(rp_p) <- c("gt", "case", "mean_fold_change", "SD" )
rp_p[ , assay := "Protein"]

rpx <- rbind(rp_r, rp_p)
rpx[ , case_assay := paste0(case, ", ", assay)]

rpx$case_assay <- factor(rpx$case_assay, levels = c("controls, RNA", "subject, RNA", 
                                                    "controls, Protein","subject, Protein" ) )

rpx$gt <- factor(rpx$gt, levels = c("nDNA-encoded", "mtDNA-encoded") )

fill_assay <- c(  "controls, RNA" = "gray80",  
                  "subject, RNA" = "#A6CEE3",
                  "controls, Protein" = "gray80", 
                  "subject, Protein" = "#FB9A99")

fig <- ggplot(rpx, aes(case_assay, mean_fold_change, fill = case_assay))+
  geom_bar(stat='identity', color = "black")+
  geom_errorbar(aes(ymin=mean_fold_change, ymax=mean_fold_change+SD), width=.2,
                position=position_dodge(.9))+
  ggtitle("Respiratoty chain complex &\nmitochondrial ribosome subunits")+
  ylab("mean fold change") +
  theme_classic()+
  scale_fill_manual( values = fill_assay)+
  facet_wrap( ~gt, ncol = 2) +

  theme(plot.title = element_text( size=12, face="bold", hjust = 0.5),
        axis.title.y= element_text( size=12, face="bold"),
        axis.title.x= element_blank(),
        axis.text.x = element_text(face="bold",  size=10),
        axis.text.y= element_text(face="bold",  size=10),
        legend.title = element_blank(),
        legend.position = "bottom",
        strip.background = element_blank(),
        strip.text = element_text(face="bold", size=12),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

fig

pdf(snakemake@output$fig, 
    width = 5, height =5,  useDingbats=FALSE )
print(fig) 
dev.off()





