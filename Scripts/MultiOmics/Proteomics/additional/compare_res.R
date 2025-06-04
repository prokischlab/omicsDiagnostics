# Investigate proteomics results with rett
source("src/functions/LIMMA/limma_functions.R")
source("src/functions/plots.R")
sa <- fread('/data/agprokisch/working/data/proteomics/raw_data/proteomics_annotation.tsv')
unique(sa$GROUP)
sa[ GROUP == "", GROUP := NA]
sa[ CATEGORY %in% c("IIb", "IIc"), SOLVED := "NO"]
sa <- sa[USE_FOR_PROTEOMICS_PAPER == T | ( USE == T & GROUP == "RETT_RETT_LIKE" ) ]
sa[ is.na(KNOWN_MUTATION), KNOWN_MUTATION := "UNSOLVED"]

saSolved <- sa[SOLVED == "YES", c("PROTEOME_ID", "KNOWN_MUTATION"  ) ]


res_unipr <- readRDS('/data/agprokisch/working/data/proteomics/processed_results/protrider/paper_results_10plex_savee.rds') %>% as.data.table()
paste("Number of samples:",  uniqueN(res_unipr$PROTEOME_ID) )
paste("Number of proteins analysed:",  uniqueN(res_unipr$geneID) ) # 7632
res_unipr[ , PROTEIN_outlier := PROTEIN_PADJ < 0.1]
os <- res_unipr[ PROTEIN_outlier == T , .N, by = PROTEOME_ID ]
paste("Median:",  median(os$N) ) # 7
plotAberrantProteinPerSample(os)
setnames(os, "N", "Outmliers_unipr")
saX <- merge(sa, os, by.x = "SAMPLE_ID", by.y = "PROTEOME_ID", all.x = T)


res_unipr <- merge( res_unipr, saSolved, by = "PROTEOME_ID", all.x = T)
res_unipr[ is.na(KNOWN_MUTATION ), KNOWN_MUTATION := "unsolved" ]
res_unipr[ , causal_gene := geneID == KNOWN_MUTATION]



res_ens <- readRDS('/data/agprokisch/working/data/proteomics/processed_results/protrider/results_10plex_savee.rds') %>% as.data.table()

paste("Number of samples:",  uniqueN(res_ens$PROTEOME_ID) )
paste("Number of proteins analysed:",  uniqueN(res_ens$geneID) ) # 8165
res_ens[ , PROTEIN_outlier := PROTEIN_PADJ < 0.1]
os <- res_ens[ PROTEIN_outlier == T , .N, by = PROTEOME_ID ]
paste("Median:",  median(os$N) ) # 8
plotAberrantProteinPerSample(os)
setnames(os, "N", "Outmliers_ens")
saX <- merge(saX, os, by = "PROTEOME_ID", all.x = T)

res_ens <- merge( res_ens, saSolved, by = "PROTEOME_ID", all.x = T)
res_ens[ is.na(KNOWN_MUTATION ), KNOWN_MUTATION := "unsolved" ]
res_ens[ , causal_gene := geneID == KNOWN_MUTATION]


res_ens_rett <- readRDS('/data/agprokisch/working/data/proteomics/processed_results/protrider/paper_rett_results_savee.rds') %>% as.data.table()

paste("Number of samples:",  uniqueN(res_ens_rett$PROTEOME_ID) )
paste("Number of proteins analysed:",  uniqueN(res_ens_rett$geneID) ) # 8042
res_ens_rett[ , PROTEIN_outlier := PROTEIN_PADJ < 0.1]
os <- res_ens_rett[ PROTEIN_outlier == T , .N, by = PROTEOME_ID ]
paste("Median:",  median(os$N) ) # 8
plotAberrantProteinPerSample(os)
setnames(os, "N", "Outmliers_ens_rett")
saX <- merge(saX, os, by = "PROTEOME_ID", all.x = T)

res_ens_rett <- merge( res_ens_rett, saSolved, by = "PROTEOME_ID", all.x = T)
res_ens_rett[ is.na(KNOWN_MUTATION ), KNOWN_MUTATION := "unsolved" ]
res_ens_rett[ , causal_gene := geneID == KNOWN_MUTATION]

# saY <- saX[ !(PROTEOME_ID %in% unique(res_ens_rett$PROTEOME_ID) )]

saX[ is.na( Outmliers_unipr), Outmliers_unipr := 0]
median(saX[PLEX == "P10"]$Outmliers_unipr)
saX[ is.na( Outmliers_ens), Outmliers_ens := 0]
median(saX[PLEX == "P10"]$Outmliers_ens)
saX[ is.na( Outmliers_ens_rett), Outmliers_ens_rett := 0]
median(saX[PLEX == "P10"]$Outmliers_ens_rett)
median(saX[AFFECTED == 0 ]$Outmliers_ens_rett)
median(saX[AFFECTED == 1]$Outmliers_ens_rett)

saX[AFFECTED == 0 , affected_status := "Unaffected"]
saX[AFFECTED == 1 , affected_status := "Affected"]

saX$affected_status <- factor(saX$affected_status, levels = c("Unaffected",  "Affected" ) )

saX[ , .N, by = affected_status]

ggplot( saX, aes(affected_status, (Outmliers_ens_rett + 0.1), color = affected_status ))+
  geom_hline(yintercept = median(saX[PLEX == "P10" ]$Outmliers_ens_rett), linetype = "dashed"  )+ 
  geom_boxplot() +
  geom_quasirandom( size = 3, alpha = 0.5) +
  #scale_color_ptol() +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", function(x) round(10^x, 0) )) +  # math_format(10^.x)
  annotation_logticks(sides = "l") +
  scale_color_colorblind()+
  ylab("# outliers per sample") +
  stat_compare_means() + 
  theme_classic()+
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.text = element_text(face = "bold", size  = 12),
        axis.title.y = element_text(face = "bold", size  = 12))



ggplot( saX[PLEX == "P11"], aes(affected_status, (Outmliers_ens_rett + 0.1), color = affected_status ))+
  geom_hline(yintercept = median(saX[PLEX == "P10" ]$Outmliers_ens_rett), linetype = "dashed"  )+ 
  geom_boxplot() +
  geom_quasirandom( size = 3, alpha = 0.5) +
  #scale_color_ptol() +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", function(x) round(10^x, 0) )) +  # math_format(10^.x)
  annotation_logticks(sides = "l") +
  scale_color_colorblind()+
  ylab("# outliers per sample") +
  stat_compare_means() + 
  theme_classic()+
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.text = element_text(face = "bold", size  = 12),
        axis.title.y = element_text(face = "bold", size  = 12))



saX <- saX[ order(Outmliers_ens_rett)] 
saX[ , SampleRank := 1 : .N]


saX[is.na(CATEGORY), CATEGORY := "unsolved"]
ggplot(saX[ !is.na(Outmliers_ens) & Outmliers_ens != 0], aes(Outmliers_ens, Outmliers_ens_rett )) +
  geom_point(aes(color = CATEGORY ), size = 3, alpha = 0.6)+
  #geom_label_repel(aes(label = geneID)) +
  geom_abline()+
  coord_fixed(xlim = c(0, 250), ylim = c(0, 250))+
  xlab("# outliers per sample\nMain cohort") +
  ylab("# outliers per sample\nExtended cohort") +
  stat_cor(method = "pearson") +
  scale_color_ptol() +
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5, size=14),
        axis.title.x= element_text( size=10, margin = NULL,face="bold"),
        axis.title.y= element_text( size=10, margin = NULL,face="bold"),
        axis.text.x = element_text(face="bold",  size=10),
        axis.text.y = element_text(face="bold",  size=10),
        legend.title = element_text(face="bold",  size=10),
        legend.position = "right")




res_ens_rett <- res_ens_rett[ order(PROTEIN_ZSCORE)]
res_ens_rett[ , SampleRank := 1 : .N, by = geneID ]

res_ens <- res_ens[ order(PROTEIN_ZSCORE)]
res_ens[ , SampleRank := 1 : .N, by = geneID ]


os <- res_ens[ PROTEIN_outlier == T , .N, by = PROTEOME_ID ]
os1 <- res_ens_rett[ PROTEOME_ID %in% unique(res_ens$PROTEOME_ID) &  PROTEIN_outlier == T , .N, by = PROTEOME_ID ]

os[ , Cohort := "Main cohort\n147 samples"]
os1[ , Cohort := "Extended cohort\n231 samples"]


os <- rbind(os, os1)
os$Cohort <- factor(os$Cohort, levels = c("Main cohort\n147 samples",  "Extended cohort\n231 samples"))
rm(os , os1)

ggplot( os, aes(Cohort, (N + 0.1), color = Cohort ))+
  geom_hline(yintercept = median(os[Cohort == "Main cohort\n147 samples" ]$N), linetype = "dashed"  )+ 
  geom_boxplot() +
  geom_quasirandom( size = 3, alpha = 0.5) +
  #scale_color_ptol() +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", function(x) round(10^x, 0) )) +  # math_format(10^.x)
  annotation_logticks(sides = "l") +
  scale_color_ptol()+
  ylab("# outliers per sample") +
  stat_compare_means() + 
  theme_classic()+
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.text = element_text(face = "bold", size  = 12),
        axis.title.y = element_text(face = "bold", size  = 12))

###################



ggplot(res_ens_rett[geneID == "MECP2" ], aes(SampleRank, PROTEIN_ZSCORE )) +
  geom_hline(yintercept = 2, linetype= "dashed", color = "grey")+
  geom_hline(yintercept = -2, linetype= "dashed", color = "grey")+
  
  geom_point(aes(color= causal_gene), size = 3, alpha = 0.6) + 
  geom_point(data = res_ens_rett[causal_gene == T & geneID == "MECP2"] , aes(SampleRank, PROTEIN_ZSCORE, color= causal_gene), size = 2.5, alpha = 0.6) + 
  
  ggtitle("MECP2") +
  xlab("Sample rank") + 
  ylab("Protein Z-score") +
  scale_color_manual( values = c("grey80", "darkred"))+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5, size=14, face="bold"),
        axis.title.x= element_text( size=10, margin = NULL,face="bold"),
        axis.title.y= element_text( size=10, margin = NULL,face="bold"),
        axis.text.x = element_text(face="bold",  size=10),
        axis.text.y = element_text(face="bold",  size=10),
        legend.title = element_text(face="bold",  size=10)) 

ggplot(res_ens_rett[geneID == "MECP2" ], aes(PROTEIN_ZSCORE )) +
  geom_histogram(aes(fill= causal_gene), alpha = 0.5, bins = 60) +
  # geom_density(aes(fill= causal_gene), alpha = 0.5) +
  scale_fill_manual( values = c("grey80", "darkred"))+
  theme_bw()


########################################

# Causal genes 

res_uniprS <- res_unipr[ causal_gene == T]
res_uniprS[ , type := "uniprot"]
res_ensS <- res_ens[ causal_gene == T]
res_ensS[ , type := "Main cohort\n147 samples"]
res_ens_rettS <- res_ens_rett[ causal_gene == T]
res_ens_rettS[ , type := "Extended cohort\n231 samples"]

uniqueN(res_ens$PROTEOME_ID)
uniqueN(res_ens_rett$PROTEOME_ID)

dfS <- rbind( res_ensS, res_ens_rettS)
dfS <- merge( dfS, sa[ , c("PROTEOME_ID", "KNOWN_MUTATION", "CATEGORY"  ) ] , by = "PROTEOME_ID", all.x = T)
dfS[ , validated := F ]
dfS[ abs(PROTEIN_ZSCORE) >= 2  | PROTEIN_PVALUE < 0.05 , validated := T ]

dfS[ , outlier_type := "non outlier"]
dfS[validated == T  , outlier_type := "nominal significance"]
dfS[PROTEIN_outlier == T  , outlier_type := "outlier"]


dfS$type <- factor(dfS$type, levels = c("Main cohort\n147 samples",  "Extended cohort\n231 samples"))

ggplot(dfS[ geneID %in% unique( res_ensS$geneID)], aes(type, PROTEIN_ZSCORE )) +
         geom_hline(yintercept = 2, linetype= "dashed", color = "grey30")+
         geom_hline(yintercept = -2, linetype= "dashed", color = "grey30")+
         geom_line(aes(group = PROTEOME_ID), color = "grey60")+
         geom_point(aes(fill= outlier_type), size = 3, shape = 21, alpha = 0.5) + 
         xlab("sample size") + 
         ylab("Protein Z-score") +
         scale_fill_manual( values = c("darkorange", "grey50", "darkred"))+
         theme_bw()+
  theme(plot.title = element_text(hjust = 0.5, size=14),
        axis.title.x= element_text( size=10, margin = NULL,face="bold"),
        axis.title.y= element_text( size=10, margin = NULL,face="bold"),
        axis.text.x = element_text(face="bold",  size=10),
        axis.text.y = element_text(face="bold",  size=10),
        legend.title = element_text(face="bold",  size=10),
        legend.position = "top")

dfS[ , .N, by = .(CATEGORY, validated )]
dfS[ , .N, by = .(CATEGORY, PROTEIN_outlier )]


df <- merge( res_ensS, res_ens_rettS, by = c("PROTEOME_ID", "geneID" ))
df <- merge( df, sa[ , c("PROTEOME_ID", "KNOWN_MUTATION", "CATEGORY"  ) ] , by = "PROTEOME_ID", all.x = T)
ggplot(df[!is.na(CATEGORY)], aes(PROTEIN_ZSCORE.x, PROTEIN_ZSCORE.y )) +
  geom_vline(xintercept = 0, color = "grey30") +
  geom_hline(yintercept = 0, color = "grey30") +
  geom_vline(xintercept = -2, linetype = "dashed", color = "grey30") +
  geom_hline(yintercept = -2, linetype = "dashed", color = "grey30") +
  geom_point(aes(color = CATEGORY ), size = 3, alpha = 0.6)+
  #geom_label_repel(aes(label = geneID)) +
  geom_abline()+
  coord_fixed(xlim = c(-13, 4), ylim = c(-13, 4))+
  xlab("Z-score Main cohort") +
  ylab("Z-score Extended cohort") +
  stat_cor(method = "pearson") +
  scale_color_ptol() +
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5, size=14),
        axis.title.x= element_text( size=10, margin = NULL,face="bold"),
        axis.title.y= element_text( size=10, margin = NULL,face="bold"),
        axis.text.x = element_text(face="bold",  size=10),
        axis.text.y = element_text(face="bold",  size=10),
        legend.title = element_text(face="bold",  size=10),
        legend.position = "top")





############################
###############################

# MECP2


prot_iterated_a <- readRDS('/data/agprokisch/working/data/proteomics/processed_results/protrider/paper_rett_iterated_a.rds') %>% as.data.table()
prot_iterated_a <- prot_iterated_a[featureID %in% unique(saSolved$KNOWN_MUTATION) ]
prot_iterated_a[ , Iteration := "a"]

prot_iterated_b <- readRDS('/data/agprokisch/working/data/proteomics/processed_results/protrider/paper_rett_iterated_b.rds') %>% as.data.table()
prot_iterated_b <- prot_iterated_b[featureID %in% unique(saSolved$KNOWN_MUTATION) ]
prot_iterated_b[ , Iteration := "b"]

prot_iterated_c <- readRDS('/data/agprokisch/working/data/proteomics/processed_results/protrider/paper_rett_iterated_c.rds') %>% as.data.table()
prot_iterated_c <- prot_iterated_c[featureID %in% unique(saSolved$KNOWN_MUTATION) ]
prot_iterated_c[ , Iteration := "c"]


prot_iterated <- rbind(prot_iterated_a, prot_iterated_b, prot_iterated_c)
prot_iterated[ , design := "ae"]
rm(prot_iterated_a, prot_iterated_b, prot_iterated_c)

setnames(prot_iterated,
         c( "featureID", "sampleID", "preprocessed_raw", "preprocessed_expected", "normalized", "zScore", "fc", "log2fc",
            "pValue", "padjust", "aberrant" ),
         c( "geneID", "PROTEOME_ID", "PROTEIN_LOG2INT_RAW", "PROTEIN_LOG2INT_EXP", "PROTEIN_LOG2INT","PROTEIN_ZSCORE", "PROTEIN_FC","PROTEIN_LOG2FC",
            "PROTEIN_PVALUE" , "PROTEIN_PADJ", "PROTEIN_outlier" ))

prot_iterated <- merge( prot_iterated, saSolved, by = "PROTEOME_ID", all.x = T)
prot_iterated[ is.na(KNOWN_MUTATION ), KNOWN_MUTATION := "unsolved" ]
prot_iterated[ , causal_gene := geneID == KNOWN_MUTATION]
prot_iterated <- prot_iterated[ order(PROTEIN_ZSCORE)]
prot_iterated[ , SampleRank := 1 : .N, by = .(IterationN, geneID, Iteration) ]
prot_iterated[ , validated := F ]
prot_iterated[ abs(PROTEIN_ZSCORE) >= 2  | PROTEIN_PVALUE < 0.05 , validated := T ]


solved_mecp2 <-  prot_iterated[ PROTEOME_ID %in% saSolved[KNOWN_MUTATION == "MECP2" ]$PROTEOME_ID ]
solved_mecp2_tp <- solved_mecp2[ causal_gene == T]

solved_mecp2_tp[validated == T, TP_v := .N, by = .(IterationN, Iteration)]
solved_mecp2_tp[ is.na(TP_v), TP_v := 0]
solved_mecp2_tp[PROTEIN_outlier == T, TP_o := .N, by = .(IterationN, Iteration)]
solved_mecp2_tp[ is.na(TP_o), TP_o := 0]

solved_mecp2_tp[  , Recall_v := TP_v / IterationN]
solved_mecp2_tp[  , Recall_o := TP_o / IterationN]

mecp2_recall <- solved_mecp2_tp[ , c( "IterationN", "Iteration", "Recall_v", "Recall_o")]
mecp2_recall <- mecp2_recall[!duplicated(mecp2_recall)]
mecp2_recall <- mecp2_recall[ order(Recall_v, Recall_o, decreasing = T )]
#mecp2_recall <- mecp2_recall[ !duplicated(IterationN)]

ggplot(mecp2_recall, aes(IterationN, Recall_v )) +
  #geom_point(aes(color = Iteration)) + 
  geom_smooth( color= "darkorange"  ) + 
  #geom_line( color = "darkorange", linewidth = 1.1) + 
  #geom_line( data = mecp2_recall, aes(IterationN, Recall_o ), color = "darkred", linewidth = 1.1) + 
  geom_smooth( data = mecp2_recall, aes(IterationN, Recall_o ), color = "darkred", linewidth = 1.1) + 
  ggtitle("MECP2 recall outrider2") +
  xlab("# samples with deffect") + 
  ylab("Recall") +
  #scale_color_manual( values = c("grey80", "darkred")) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, size=14),
        axis.title.x= element_text( size=10, margin = NULL,face="bold"),
        axis.title.y= element_text( size=10, margin = NULL,face="bold"),
        axis.text.x = element_text(face="bold",  size=10),
        axis.text.y = element_text(face="bold",  size=10),
        legend.title = element_blank()  )

mecp2_recall_ae <- mecp2_recall

ggplot(mecp2_recall_ae[Recall_v != 0], aes(IterationN, Recall_v )) +
  #geom_point(aes(color = Iteration)) + 
  geom_smooth( color= "darkorange"  ) + 
  geom_smooth( data = mecp2_recall_ae[Recall_o != 0], aes(IterationN, Recall_o ), color = "darkred", linewidth = 1.1) + 
  ggtitle("MECP2 recall outrider2") +
  xlab("# samples with deffect") + 
  ylab("Recall") +
  #scale_color_manual( values = c("grey80", "darkred")) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, size=14),
        axis.title.x= element_text( size=10, margin = NULL,face="bold"),
        axis.title.y= element_text( size=10, margin = NULL,face="bold"),
        axis.text.x = element_text(face="bold",  size=10),
        axis.text.y = element_text(face="bold",  size=10),
        legend.title = element_blank()  )



mecp2 <- prot_iterated[geneID == "MECP2" ]

vlidated_cases <- unique(mecp2[abs(PROTEIN_ZSCORE) >= 2  | PROTEIN_PVALUE < 0.05 ]$PROTEOME_ID)




ggplot(mecp2, aes(SampleRank, PROTEIN_ZSCORE )) +
  geom_hline(yintercept = 2, linetype= "dashed", color = "grey")+
  geom_hline(yintercept = -2, linetype= "dashed", color = "grey")+
  geom_point(aes(color= causal_gene)) + 
  ggtitle("MECP2") +
  xlab("Sample rank") + 
  ylab("Protein Z-score") +
  scale_color_manual( values = c("grey80", "darkred"))+
  
  # geom_label_repel(data= samp[outlier_class != "non_outlier" & gene_class == "rare pot. biallelic variants"], 
  #                  mapping=aes(rank_protein, PROTEIN_INT, label= SAMPLE_ID ), 
  #                  box.padding = unit(0.3, "lines"), colour= "black", 
  #                  size= 4, show.legend = F )+ 
  # 
  
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5, size=14),
        axis.title.x= element_text( size=10, margin = NULL,face="bold"),
        axis.title.y= element_text( size=10, margin = NULL,face="bold"),
        axis.text.x = element_text(face="bold",  size=10),
        axis.text.y = element_text(face="bold",  size=10),
        legend.title = element_blank(),
        plot.margin = margin(0, 0, 0, 0, "cm") ) +
  facet_wrap(~IterationN)





ggplot(mecp2[ causal_gene == T & PROTEOME_ID %in%  vlidated_cases], aes(IterationN, PROTEIN_ZSCORE )) +
  geom_hline(yintercept = 2, linetype= "dashed", color = "grey30")+
  geom_hline(yintercept = -2, linetype= "dashed", color = "grey30")+
  geom_line(aes(group = PROTEOME_ID), color = "grey60")+
  geom_point(aes(color= validated)) + 
  ggtitle("MECP2") +
  xlab("# samples with deffect") + 
  ylab("Protein Z-score") +
  scale_color_manual( values = c("grey50", "darkred"))+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5, size=14),
        axis.title.x= element_text( size=10, margin = NULL,face="bold"),
        axis.title.y= element_text( size=10, margin = NULL,face="bold"),
        axis.text.x = element_text(face="bold",  size=10),
        axis.text.y = element_text(face="bold",  size=10),
        legend.title = element_text(face="bold",  size=10) )


ggplot(mecp2[ causal_gene == T & PROTEOME_ID %in%  vlidated_cases], aes(IterationN, -log10(PROTEIN_PVALUE) )) +
  geom_hline(yintercept = -log10(0.05), linetype= "dashed", color = "grey")+
  geom_line(aes(group = PROTEOME_ID), color = "grey60")+
  geom_point(aes(color= validated)) + 
  ggtitle("MECP2") +
  xlab("# samples with deffect") + 
  ylab("-log10( p-value)") +
  scale_color_manual( values = c("grey50", "darkred"))+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5, size=14),
        axis.title.x= element_text( size=10, margin = NULL,face="bold"),
        axis.title.y= element_text( size=10, margin = NULL,face="bold"),
        axis.text.x = element_text(face="bold",  size=10),
        axis.text.y = element_text(face="bold",  size=10),
        legend.title = element_text(face="bold",  size=10) )





ggplot(mecp2[ causal_gene == T & PROTEOME_ID %in%  vlidated_cases], aes(IterationN, PROTEIN_LOG2INT )) +
 # geom_hline(yintercept = -log10(0.05), linetype= "dashed", color = "grey")+
  geom_line(aes(group = PROTEOME_ID), color = "grey60")+
  geom_point(aes(color= validated)) + 
  ggtitle("MECP2") +
  xlab("# samples with deffect") + 
  ylab("Log2 Intensity") +
  scale_color_manual( values = c("grey50", "darkred"))+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5, size=14),
        axis.title.x= element_text( size=10, margin = NULL,face="bold"),
        axis.title.y= element_text( size=10, margin = NULL,face="bold"),
        axis.text.x = element_text(face="bold",  size=10),
        axis.text.y = element_text(face="bold",  size=10),
        legend.title = element_text(face="bold",  size=10) )





####################################################
prot_iterations_mecp2 <- readRDS('/data/agprokisch/working/data/proteomics/processed_results/protrider/iterate_full_dataset_MECP2.rds') %>% as.data.table()
prot_iterations_mecp2 <- merge( prot_iterations_mecp2, sa[ , c("PROTEOME_ID", "KNOWN_MUTATION", "CATEGORY"  ) ] , by = "PROTEOME_ID", all.x = T)
prot_iterations_mecp2[ is.na(KNOWN_MUTATION ), KNOWN_MUTATION := "unsolved" ]
prot_iterations_mecp2[ , causal_gene := geneID == KNOWN_MUTATION]
prot_iterations_mecp2[ , validated := F ]
prot_iterations_mecp2[ abs(PROTEIN_ZSCORE) >= 2  | PROTEIN_PVALUE < 0.05 , validated := T ]




ggplot(prot_iterations_mecp2[ causal_gene == T], aes(as.factor(SampleSize),  -log10(PROTEIN_PVALUE) )  ) +
  geom_hline(yintercept = -log10(0.05), linetype= "dashed", color = "grey")+
  # geom_line(aes(group = PROTEOME_ID), color = "grey60")+
  geom_boxplot(aes(color= causal_gene)) + 
  ggtitle("MECP2") +
  xlab("# sample size") + 
  ylab("-log10( p-value)") +
  scale_color_manual( values = c("grey50", "darkred"))+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5, size=14),
        axis.title.x= element_text( size=10, margin = NULL,face="bold"),
        axis.title.y= element_text( size=10, margin = NULL,face="bold"),
        axis.text.x = element_text(face="bold",  size=10),
        axis.text.y = element_text(face="bold",  size=10),
        legend.title = element_text(face="bold",  size=10) )



ggplot(prot_iterations_mecp2[ causal_gene == T], aes(as.factor(SampleSize),  -log10(PROTEIN_PVALUE))  ) +
  geom_hline(yintercept = -log10(0.05), linetype= "dashed", color = "grey")+
  # geom_line(aes(group = PROTEOME_ID), color = "grey60")+
  geom_boxplot(aes(color= causal_gene)) + 
  ggtitle("MECP2") +
  xlab("# sample size") + 
  ylab("-log10( p-value)") +
  scale_color_manual( values = c("grey50", "darkred"))+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5, size=14),
        axis.title.x= element_text( size=10, margin = NULL,face="bold"),
        axis.title.y= element_text( size=10, margin = NULL,face="bold"),
        axis.text.x = element_text(face="bold",  size=10),
        axis.text.y = element_text(face="bold",  size=10),
        legend.title = element_text(face="bold",  size=10) )




ggplot(prot_iterations_mecp2[ PROTEOME_ID %in% unique(prot_iterations[validated == T ]$PROTEOME_ID) ], aes(as.factor(SampleSize),   -log10(PROTEIN_PVALUE) )  ) +
  geom_boxplot(aes(color= causal_gene)) + 
  geom_hline(yintercept = -log10(0.05), linetype= "dashed", color = "grey40")+
  ggtitle("MECP2") +
  xlab("# sample size") + 
  ylab("-log10( p-value)") +
  scale_color_manual( values = c("grey50", "darkred"))+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5, size=14),
        axis.title.x= element_text( size=10, margin = NULL,face="bold"),
        axis.title.y= element_text( size=10, margin = NULL,face="bold"),
        axis.text.x = element_text(face="bold",  size=10),
        axis.text.y = element_text(face="bold",  size=10),
        legend.title = element_text(face="bold",  size=10) )+
  facet_wrap(~causal_gene, nrow = 2)




ggplot(prot_iterations_mecp2[ PROTEOME_ID %in% unique(prot_iterations_mecp2[validated == T ]$PROTEOME_ID) ], aes((SampleSize),   -log10(PROTEIN_PVALUE) )  ) +
  geom_smooth(aes(color= causal_gene)) + 
  geom_hline(yintercept = -log10(0.05), linetype= "dashed", color = "grey40")+
  ggtitle("MECP2") +
  xlab("# sample size") + 
  ylab("-log10( p-value)") +
  scale_color_manual( values = c("grey50", "darkred"))+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5, size=14),
        axis.title.x= element_text( size=10, margin = NULL,face="bold"),
        axis.title.y= element_text( size=10, margin = NULL,face="bold"),
        axis.text.x = element_text(face="bold",  size=10),
        axis.text.y = element_text(face="bold",  size=10),
        legend.title = element_text(face="bold",  size=10) ) # + facet_wrap(~causal_gene, nrow = 2)


prot_iterations_mecp2 <- readRDS('/data/agprokisch/working/data/proteomics/processed_results/protrider/iterate_full_dataset_MECP2.rds') %>% as.data.table()
prot_iterations_mecp2 <- merge( prot_iterations_mecp2, sa[ , c("PROTEOME_ID", "KNOWN_MUTATION", "CATEGORY"  ) ] , by = "PROTEOME_ID", all.x = T)
prot_iterations_mecp2[ is.na(KNOWN_MUTATION ), KNOWN_MUTATION := "unsolved" ]
prot_iterations_mecp2[ , causal_gene := geneID == KNOWN_MUTATION]
prot_iterations_mecp2[ , validated := F ]
prot_iterations_mecp2[ abs(PROTEIN_ZSCORE) >= 2  | PROTEIN_PVALUE < 0.05 , validated := T ]




solved_mecp2 <-  prot_iterations_mecp2[ PROTEOME_ID %in% saSolved[KNOWN_MUTATION == "MECP2" ]$PROTEOME_ID ]
solved_mecp2_tp <- solved_mecp2[ causal_gene == T]

solved_mecp2_tp[ , TP := uniqueN(PROTEOME_ID), by = .(SampleSize, Iteration)]

solved_mecp2_tp[validated == T, TP_v := .N, by = .(SampleSize, Iteration)]
solved_mecp2_tp[ is.na(TP_v), TP_v := 0]

solved_mecp2_tp[  , Recall := TP_v / TP]


mecp2_recall <- solved_mecp2_tp[ , c( "SampleSize" ,"Iteration", "Recall")]
mecp2_recall <- mecp2_recall[!duplicated(mecp2_recall)]
mecp2_recall <- mecp2_recall[ order(Recall,  decreasing = T )]
mecp2_recall <- mecp2_recall[ !duplicated(Iteration, SampleSize)]

mecp2_recall <- mecp2_recall[ order(SampleSize, Iteration )]


ggplot(mecp2_recall[Recall != 0 ], aes(SampleSize, Recall )) +
  #geom_point( color = "darkorange" ) + 
  geom_smooth( color= "darkorange"  ) + 
  ggtitle("MECP2 recall") +
  xlab("sample size") + 
  ylab("Recall") +
  scale_color_manual( values = c("darkorange" )) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, size=14),
        axis.title.x= element_text( size=10, margin = NULL,face="bold"),
        axis.title.y= element_text( size=10, margin = NULL,face="bold"),
        axis.text.x = element_text(face="bold",  size=10),
        axis.text.y = element_text(face="bold",  size=10),
        legend.title = element_blank()  )


ggplot(mecp2_recall[Recall != 0 ], aes(SampleSize, Recall )) +
  geom_point( color = "darkorange" ) + 
  # geom_smooth( color= "darkorange"  ) + 
  ggtitle("MECP2 recall") +
  xlab("sample size") + 
  ylab("Recall") +
  scale_color_manual( values = c("darkorange" )) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, size=14),
        axis.title.x= element_text( size=10, margin = NULL,face="bold"),
        axis.title.y= element_text( size=10, margin = NULL,face="bold"),
        axis.text.x = element_text(face="bold",  size=10),
        axis.text.y = element_text(face="bold",  size=10),
        legend.title = element_blank()  )



ggplot(mecp2_recall, aes(as.factor(SampleSize), Recall )) +
  geom_boxplot( color = "darkorange" ) + 
  # geom_smooth( color= "darkorange"  ) + 
  ggtitle("MECP2 recall") +
  xlab("sample size") + 
  ylab("Recall") +
  scale_color_manual( values = c("darkorange" )) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, size=14),
        axis.title.x= element_text( size=10, margin = NULL,face="bold"),
        axis.title.y= element_text( size=10, margin = NULL,face="bold"),
        axis.text.x = element_text(face="bold",  size=10),
        axis.text.y = element_text(face="bold",  size=10),
        legend.title = element_blank()  )



##################################################


prot_iterations <- readRDS('/data/agprokisch/working/data/proteomics/processed_results/protrider/iterate_full_dataset_causal.rds') %>% as.data.table()
prot_iterations <- merge( prot_iterations, sa[ , c("PROTEOME_ID", "KNOWN_MUTATION", "CATEGORY"  ) ] , by = "PROTEOME_ID", all.x = T)
prot_iterations[ is.na(KNOWN_MUTATION ), KNOWN_MUTATION := "unsolved" ]
prot_iterations[ , causal_gene := geneID == KNOWN_MUTATION]
prot_iterations[ , validated := F ]
prot_iterations[ abs(PROTEIN_ZSCORE) >= 2  | PROTEIN_PVALUE < 0.05 , validated := T ]


ggplot(prot_iterations[ PROTEOME_ID %in% unique(prot_iterations[validated == T ]$PROTEOME_ID) ], aes((SampleSize),   -log10(PROTEIN_PVALUE) )  ) +
  geom_smooth(aes(color= causal_gene)) + 
  geom_hline(yintercept = -log10(0.05), linetype= "dashed", color = "grey40")+
  ggtitle("Causal genes") +
  xlab("# sample size") + 
  ylab("-log10( p-value)") +
  scale_color_manual( values = c("grey50", "darkred"))+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5, size=14),
        axis.title.x= element_text( size=10, margin = NULL,face="bold"),
        axis.title.y= element_text( size=10, margin = NULL,face="bold"),
        axis.text.x = element_text(face="bold",  size=10),
        axis.text.y = element_text(face="bold",  size=10),
        legend.title = element_text(face="bold",  size=10) ) # + facet_wrap(~causal_gene, nrow = 2)







solved <-  prot_iterations[ PROTEOME_ID %in% saSolved$PROTEOME_ID ]
solved <- solved[ causal_gene == T]

solved[ , TP := uniqueN(PROTEOME_ID), by = .(SampleSize, Iteration)]

solved[validated == T, TP_v := .N, by = .(SampleSize, Iteration)]
solved[ is.na(TP_v), TP_v := 0]
solved[  , Recall := TP_v / TP]


solved_recall <- solved[ , c( "SampleSize" ,"Iteration", "Recall")]
solved_recall <- solved_recall[!duplicated(solved_recall)]
solved_recall <- solved_recall[ order(Recall,  decreasing = T )]
solved_recall <- solved_recall[ !duplicated(SampleSize, Iteration)]

solved_recall <- solved_recall[ order(SampleSize, Iteration )]


ggplot(solved_recall[Recall != 0 ], aes(SampleSize, Recall )) +
  #geom_point( color = "darkorange" ) + 
  geom_smooth( color= "darkred"  ) + 
  ggtitle("Causal gene recall") +
  xlab("sample size") + 
  ylab("Recall") +
  scale_color_manual( values = c("darkorange" )) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, size=14),
        axis.title.x= element_text( size=10, margin = NULL,face="bold"),
        axis.title.y= element_text( size=10, margin = NULL,face="bold"),
        axis.text.x = element_text(face="bold",  size=10),
        axis.text.y = element_text(face="bold",  size=10),
        legend.title = element_blank()  )



ggplot(solved_recall[Recall != 0 ], aes(SampleSize, Recall )) +
  geom_point( color = "darkorange" ) + 
  # geom_smooth( color= "darkorange"  ) + 
  ggtitle("Causal gene recall") +
  xlab("sample size") + 
  ylab("Recall") +
  scale_color_manual( values = c("darkorange" )) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, size=14),
        axis.title.x= element_text( size=10, margin = NULL,face="bold"),
        axis.title.y= element_text( size=10, margin = NULL,face="bold"),
        axis.text.x = element_text(face="bold",  size=10),
        axis.text.y = element_text(face="bold",  size=10),
        legend.title = element_blank()  )



ggplot(solved_recall, aes(as.factor(SampleSize), Recall )) +
  geom_boxplot( color = "darkorange" ) + 
  # geom_smooth( color= "darkorange"  ) + 
  ggtitle("Causal gene recall") +
  xlab("sample size") + 
  ylab("Recall") +
  scale_color_manual( values = c("darkorange" )) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, size=14),
        axis.title.x= element_text( size=10, margin = NULL,face="bold"),
        axis.title.y= element_text( size=10, margin = NULL,face="bold"),
        axis.text.x = element_text(face="bold",  size=10),
        axis.text.y = element_text(face="bold",  size=10),
        legend.title = element_blank()  )






####################################


prot_iterations_mecp2 <- readRDS('/data/agprokisch/working/data/proteomics/processed_results/protrider/iterate_MECP2_causal.rds') %>% as.data.table()
prot_iterations_mecp2 <- merge( prot_iterations_mecp2, sa[ , c("PROTEOME_ID", "KNOWN_MUTATION", "CATEGORY"  ) ] , by = "PROTEOME_ID", all.x = T)
prot_iterations_mecp2[ is.na(KNOWN_MUTATION ), KNOWN_MUTATION := "unsolved" ]
prot_iterations_mecp2[ , causal_gene := geneID == KNOWN_MUTATION]
prot_iterations_mecp2[ , validated := F ]
prot_iterations_mecp2[ abs(PROTEIN_ZSCORE) >= 2  | PROTEIN_PVALUE < 0.05 , validated := T ]




solved_mecp2 <-  prot_iterations_mecp2[ PROTEOME_ID %in% saSolved[KNOWN_MUTATION == "MECP2" ]$PROTEOME_ID ]
solved_mecp2_tp <- solved_mecp2[ causal_gene == T]

solved_mecp2_tp[ , TP := uniqueN(PROTEOME_ID), by = .(SampleSize, Iteration)]

solved_mecp2_tp[validated == T, TP_v := .N, by = .(SampleSize, Iteration)]
solved_mecp2_tp[ is.na(TP_v), TP_v := 0]

solved_mecp2_tp[  , Recall := TP_v / TP]


mecp2_recall <- solved_mecp2_tp[ , c( "SampleSize" ,"Iteration", "Recall")]
mecp2_recall <- mecp2_recall[!duplicated(mecp2_recall)]
mecp2_recall <- mecp2_recall[ order(Recall,  decreasing = T )]
mecp2_recall <- mecp2_recall[ !duplicated(Iteration, SampleSize)]

mecp2_recall <- mecp2_recall[ order(SampleSize, Iteration )]


ggplot(mecp2_recall[Recall != 0 ], aes(SampleSize, Recall )) +
  #geom_point( color = "darkorange" ) + 
  geom_smooth( color= "darkorange"  ) + 
  ggtitle("MECP2 recall") +
  xlab("sample size") + 
  ylab("Recall") +
  scale_color_manual( values = c("darkorange" )) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, size=14),
        axis.title.x= element_text( size=10, margin = NULL,face="bold"),
        axis.title.y= element_text( size=10, margin = NULL,face="bold"),
        axis.text.x = element_text(face="bold",  size=10),
        axis.text.y = element_text(face="bold",  size=10),
        legend.title = element_blank()  )



colnames(mecp2_recall_ae)[1] <- "SampleSize"
colnames(mecp2_recall_ae)[3] <- "Recall"
mecp2_recall_ae[ , type := "AE + outlier test"]
mecp2_recall[ , Recall_o := 0]
mecp2_recall[ , type := "Outlier test"]


mecp2_recallX <- rbind(mecp2_recall_ae, mecp2_recall)


ggplot(mecp2_recallX[Recall != 0 ], aes(SampleSize, Recall )) +
  #geom_point( color = "darkorange" ) + 
  geom_smooth( aes(color = type) ) + 
  ggtitle("MECP2 recall") +
  xlab("sample size") + 
  ylab("Recall") +
  #scale_color_manual( values = c("darkorange" )) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, size=14),
        axis.title.x= element_text( size=10, margin = NULL,face="bold"),
        axis.title.y= element_text( size=10, margin = NULL,face="bold"),
        axis.text.x = element_text(face="bold",  size=10),
        axis.text.y = element_text(face="bold",  size=10),
        legend.title = element_blank()  )



