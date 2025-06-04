source("kaisa_prep/Scripts/config.R")


df <- readRDS('kaisa_prep/processed_data/limma_results.rds')
sa <- fread('kaisa_prep/raw_data/proteomics_annotation.tsv')
sa <- sa[PROTEOME_ID %in% unique(df$PROTEOME_ID)]
sac <- sa[ DISEASE == "Control"]


colnames(df)
df[ , MEAN_INT := mean(PROTEIN_INT, na.rm = T ), by = "geneID"]
df[ , SD := sd(PROTEIN_INT, na.rm = T ), by = "geneID"]

df[PROTEOME_ID %in% sac$PROTEOME_ID, MEAN_INT_c := mean(PROTEIN_INT, na.rm = T ), by = "geneID"]
df[PROTEOME_ID %in% sac$PROTEOME_ID , SD_c := sd(PROTEIN_INT, na.rm = T ), by = "geneID"]

dfF <- df[!is.na(MEAN_INT_c) ]
dfF <- dfF[, c( "geneID", "MEAN_INT", "SD", "MEAN_INT_c", "SD_c") ]
dfF <- dfF[!duplicated(dfF)]

dfF[ , delta_Mean := MEAN_INT - MEAN_INT_c ]

dfF[ , perc_change :=  100 * abs(delta_Mean)  / mean(MEAN_INT, MEAN_INT_c), by = geneID ]


ggplot(dfF, aes(MEAN_INT_c, MEAN_INT))+
  geom_abline()+
  geom_point()+
  ylab("Mean abundance overall, ng") +
  xlab("Mean abundance in controls, ng") +
  geom_label_repel(data = dfF[ perc_change > 40], aes(MEAN_INT_c, MEAN_INT, label = geneID))+ 
  scale_y_log10( ) +
  scale_x_log10( ) +
  annotation_logticks(sides="bl")+ 
  coord_fixed(xlim = c(min(dfF$MEAN_INT_c, dfF$MEAN_INT), max(dfF$MEAN_INT_c, dfF$MEAN_INT)), 
              ylim = c(min(dfF$MEAN_INT_c, dfF$MEAN_INT), max(dfF$MEAN_INT_c, dfF$MEAN_INT)))+
  theme_classic()




dfX <- dfF
dfX <- dfX[!duplicated(dfX)]
dfX <- dfX[ order(MEAN_INT + SD)]
dfX <- dfX[ order(MEAN_INT )]

dfX[ , YMIN := MEAN_INT - SD ]
dfX[ , YMAX := MEAN_INT + SD ]
dfX$GENE <- dfX$geneID
dfX$GENE <- factor(dfX$GENE, levels = unique(dfX$GENE))
dfX[ , Rank := 1 : .N ]


dfX <- dfX[ order(MEAN_INT_c + SD_c)]
dfX <- dfX[ order(MEAN_INT_c )]
dfX[ , YMAX_c := MEAN_INT_c + SD_c ]

# Correct classification of amounts with a proper conditional structure
dfX[ YMAX < 0.01 , ammount_class := "picograms"]
dfX[ YMAX >= 0.01 & YMAX < 0.1 , ammount_class := "tens of picograms"]
dfX[ YMAX >= 0.1 & YMAX < 1 , ammount_class := "hundreds of picograms"]
dfX[ YMAX >= 1 & YMAX < 10 , ammount_class := "nanograms"]
dfX[ YMAX >= 10 & YMAX < 100 , ammount_class := "tens of nanograms"]
dfX[ YMAX >= 100 & YMAX < 1000 , ammount_class := "hundreds of nanograms"]
dfX[ YMAX >= 1000 , ammount_class := "micrograms"]


# Correct classification of amounts with a proper conditional structure
dfX[ YMAX_c < 0.01 , ammount_class_c := "picograms"]
dfX[ YMAX_c >= 0.01 & YMAX_c < 0.1 , ammount_class_c := "tens of picograms"]
dfX[ YMAX_c >= 0.1 & YMAX_c < 1 , ammount_class_c := "hundreds of picograms"]
dfX[ YMAX_c >= 1 & YMAX_c < 10 , ammount_class_c := "nanograms"]
dfX[ YMAX_c >= 10 & YMAX_c < 100 , ammount_class_c := "tens of nanograms"]
dfX[ YMAX_c >= 100 & YMAX_c < 1000 , ammount_class_c := "hundreds of nanograms"]
dfX[ YMAX_c >= 1000 , ammount_class_c := "micrograms"]
dfX <- dfX[!is.na(ammount_class_c)]

# Ensure that the factors are ordered correctly
dfX$ammount_class <- factor(dfX$ammount_class, 
                            levels = c("picograms", "tens of picograms", "hundreds of picograms", 
                                       "nanograms", "tens of nanograms", 
                                       "hundreds of nanograms", "micrograms"))

dfX$ammount_class_c <- factor(dfX$ammount_class_c, 
                            levels = c("picograms", "tens of picograms", "hundreds of picograms", 
                                       "nanograms", "tens of nanograms", 
                                       "hundreds of nanograms", "micrograms"))


dfZ <- dfX[ , .N, by = ammount_class]
ggplot(dfZ , aes(ammount_class, N))+
  geom_col()+
  theme_classic()+
  geom_label(aes(label = N)) +
  #annotation_logticks(sides="b")+ 
  coord_flip() +
  xlab("Ammount ranges")+ 
  ylab("Number of proteins")+ 
  theme(#axis.ticks.x = element_blank(), 
    axis.text.y = element_text(face = "bold"),
    axis.text.x = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold"),
    axis.title.x = element_text(face = "bold"))


ggplot(dfX[ , .N, by = ammount_class_c] , aes(ammount_class_c, N))+
  geom_col()+
  theme_classic()+
  geom_label(aes(label = N)) +
  #annotation_logticks(sides="b")+ 
  coord_flip() +
  xlab("Ammount ranges")+ 
  ylab("Number of proteins")+ 
  theme(#axis.ticks.x = element_blank(), 
    axis.text.y = element_text(face = "bold"),
    axis.text.x = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold"),
    axis.title.x = element_text(face = "bold"))



dfX[  YMAX < 1 , ammount_classX := "picograms"] 
dfX[  YMAX < 100 , ammount_classX := "nanograms"]
dfX[ YMAX >= 100 & YMAX < 1000 , ammount_classX := "hundreds of nanograms"]
dfX[ YMAX >= 1000 , ammount_classX := "micrograms"]

# Ensure that the factors are ordered correctly
dfX$ammount_classX <- factor(dfX$ammount_classX, 
                            levels = c("picograms",  
                                       "nanograms", 
                                       "hundreds of nanograms", "micrograms"))


#total_protein_mass <- 15000 # 15 µg , 15000 is in nanograms (ng)


ggplot(dfX[ammount_classX == "hundreds of nanograms"] , aes(GENE, MEAN_INT))+
  geom_point()+
  #scale_y_log10( ) +
  geom_errorbar(aes(ymin=MEAN_INT-SD, ymax=MEAN_INT+SD), width=.2)+
  theme_classic()+
  #annotation_logticks(sides="b")+ 
  coord_flip() +
  xlab("Proteins")+ 
  ylab("Estimated mass in nanograms")+ 
  theme(#axis.ticks.x = element_blank(), 
    #axis.text.x = element_text(angle = 90),
    axis.title.y = element_text(face = "bold"),
    axis.title.x = element_text(face = "bold"))


ggplot(dfX[ammount_class_c == "hundreds of nanograms"] , aes(GENE, MEAN_INT_c))+
  geom_point()+
  #scale_y_log10( ) +
  geom_errorbar(aes(ymin=MEAN_INT_c-SD_c, ymax=MEAN_INT_c+SD_c), width=.2)+
  theme_classic()+
  #annotation_logticks(sides="b")+ 
  coord_flip() +
  xlab("Proteins")+ 
  ylab("Estimated mass in nanograms")+ 
  theme(#axis.ticks.x = element_blank(), 
    #axis.text.x = element_text(angle = 90),
    axis.title.y = element_text(face = "bold"),
    axis.title.x = element_text(face = "bold"))





ggplot(dfX , aes(Rank, MEAN_INT))+
  geom_point()+
  scale_y_log10( ) +
  geom_errorbar(aes(ymin=MEAN_INT-SD, ymax=MEAN_INT+SD), width=.2)+
  theme_classic()+
  annotation_logticks(sides="l")+ 
  #coord_flip() +
  xlab("Protein rank")+ 
  ylab("Estimated mass in nanograms")+ 
  theme(#axis.ticks.x = element_blank(), 
        #axis.text.x = element_blank(),
        axis.title.y = element_text(face = "bold"),
        axis.title.x = element_text(face = "bold"))




ggplot(dfX, aes(Rank, MEAN_INT))+
  geom_point()+
  scale_y_log10( ) +
  scale_x_log10( ) +
  annotation_logticks(sides="bl")+ 
  geom_errorbar(aes(ymin=MEAN_INT-SD, ymax=MEAN_INT+SD), width=.2)+
  theme_classic()+

  #coord_flip() +
  xlab("Protein rank")+ 
  ylab("Estimated mass in nanograms")+ 
  theme(#axis.ticks.x = element_blank(), 
    #axis.text.x = element_blank(),
    axis.title.y = element_text(face = "bold"),
    axis.title.x = element_text(face = "bold"))


ggplot(dfX[YMIN > -5 ], aes(GENE, MEAN_INT))+
  geom_point()+
  # geom_bar(position=position_dodge(), stat="identity", 
  #          colour='black') + 
  geom_errorbar(aes(ymin=MEAN_INT-SD, ymax=MEAN_INT+SD), width=.2)+
  theme_classic()+
  theme(axis.ticks.x = element_blank())+
  # scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
  #               labels = trans_format("log10", math_format(10^.x))) +
  scale_y_log10( ) +
  coord_flip()+
  annotation_logticks(sides="b")+ 
  facet_wrap(~ammount_class, scales = "free")


saS <- sa[ , c("PROTEOME_ID", "KNOWN_MUTATION")]
saS <- saS[!duplicated(saS)]

dfS <- merge(df, saS, by = "PROTEOME_ID")
dfS[ , causal_gene := geneID == KNOWN_MUTATION] 
dfS[ is.na(causal_gene) , causal_gene := F]
colnames(df)
P1 <- ggplot(dfS[geneID %in% c("NBAS"  ) & PROTEIN_INT < 0.23 ] , aes(geneID, PROTEIN_INT))+
  geom_quasirandom(aes(color = causal_gene))+
  geom_point(data = dfS[causal_gene == T & geneID %in% c("NBAS" ) & PROTEIN_INT < 0.23 ] , aes(geneID, PROTEIN_INT, color = causal_gene), size = 3)+
  geom_point(data = dfS[ geneID %in% c(  "NBAS") & PROTEIN_INT < 0.23 ] , aes(geneID, MEAN_INT) , size = 2)+
  geom_errorbar(aes(ymin=MEAN_INT-2*SD, ymax=MEAN_INT+ 2*SD), width=.2)+
  theme_classic()+
  #scale_y_log10( ) +
  scale_color_manual(values = c( "grey60", "darkred")) +
 # annotation_logticks(sides="l")+ 
  #coord_flip() +
  xlab("")+ 
  ylab("Estimated mass (ng)")+ 
  theme(#axis.ticks.x = element_blank(), 
    #axis.text.x = element_blank(),
    axis.title.y = element_text(face = "bold"),
    axis.title.x = element_text(face = "bold"), 
    legend.position = "none")

P2 <- ggplot(dfS[geneID %in% c(  "EPG5" ) ] , aes(geneID, PROTEIN_INT))+
  geom_quasirandom(aes(color = causal_gene))+
  geom_point(data = dfS[causal_gene == T & geneID %in% c(  "EPG5") ] , aes(geneID, PROTEIN_INT, color = causal_gene), size = 3)+
  
  geom_point(data = dfS[ geneID %in% c(  "EPG5") ] , aes(geneID, MEAN_INT) , size = 2)+
  #scale_y_log10( ) +
  geom_errorbar(aes(ymin=MEAN_INT-2*SD, ymax=MEAN_INT+ 2*SD), width=.2)+
  theme_classic()+
  scale_color_manual(values = c( "grey60", "darkred")) +
  #annotation_logticks(sides="l")+ 
  #coord_flip() +
  xlab("")+ 
  ylab("Estimated mass (ng)")+ 
  theme(#axis.ticks.x = element_blank(), 
    #axis.text.x = element_blank(),
    axis.title.y = element_text(face = "bold"),
    axis.title.x = element_text(face = "bold"))

P1 | P2


