#'---
#' title: Supplementary Figure 2e, proportion of outliers with pb variants and SemSim >2 
#' author: Dmitrii Smirnov
#' wb:
#'  input:
#'  - enrichments_proportions: '`sm config["RAW_DATA"] + "/enrichment_proportions_variants.tsv"`'
#'  output:
#'  - fig: '`sm config["FIGURE_DIR"] + "/Supplementary_figures/S_Fig2_e.pdf"`'
#' output: 
#'   html_document:
#'    code_folding: hide
#'    code_download: TRUE
#'---

source('src/config.R')


# Read combined file produced on full variant table
px <- fread(snakemake@input$enrichments_proportions)
px <- px[up_down_outlier %in% c("RNA_underexpression", "Protein_underexpression", "RNA_Protein_underexpression", "non_outlier")]
px <- px[, c("up_down_outlier", "var_type", "prop", "total", "type")]


px[var_type == "potential biallelic SemSim > 2", var_type := "potential biallelic + Semantic sim > 2"] 
px[, var_type := factor(var_type, levels = c("no rare variant", "rare", 
                                             "potential_biallelic", "potential biallelic + Semantic sim > 2", 
                                             "non_coding", "synonymous", "coding", "frameshift", "splice", "stop"))]


# unique(px$up_down_outlier)
px[up_down_outlier == "RNA_Protein_underexpression", up_down_outlier :=  "RNA-and-protein"]
px[up_down_outlier == "Protein_underexpression" , up_down_outlier := "protein-only" ]
px[up_down_outlier == "RNA_underexpression" , up_down_outlier := "RNA-only" ]
px[up_down_outlier == "non_outlier", up_down_outlier :=  "non outlier"]

px[, outlier_class_label := paste0(up_down_outlier, "\n", "(n = ", total, ")") ]

# unique(px$outlier_class_label)
px$outlier_class_label <- factor(px$outlier_class_label , 
                                 levels = c( "non outlier\n(n = 217894)",
                                             "RNA-only\n(n = 180)",
                                             "protein-only\n(n = 924)", 
                                             "RNA-and-protein\n(n = 147)" ))



px <- px[ !(var_type %in% c("no rare variant", "rare") )]

#px[ type == "var_type", widt := 1]
#px[ type == "pb", widt := 0.5]

s_fig <- ggplot(px[ type != "var_type"], aes(outlier_class_label, prop)) +
  geom_bar(stat= 'identity', aes(fill = var_type)) +
  geom_hline(yintercept = 0, colour = "black") + 
  scale_fill_manual(values = variant_colors ) +
  coord_flip(ylim = c(0, 0.45)) +
  scale_y_continuous(breaks=seq(0, 0.45, 0.1), labels=scales::percent) +
  labs( y = "Proportion of underexpression outliers with rare variants")+
  # theme_classic()+
  theme(legend.position="top",  
        legend.direction = "horizontal",
        legend.title = element_blank(),
        
        
        axis.title.y = element_blank() ,
        axis.text.y = element_text(face="bold", size=12), 
        axis.title.x = element_text(face="bold", size=12) ,
        
        
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.line.x.bottom =  element_line(colour = "black"),
        plot.background = element_rect( fill = "white")) +
  guides(fill = guide_legend(nrow = 1))

#+ fig.width=7.5, fig.height=3
s_fig

pdf(snakemake@output$fig, 
    width = 7.5, height = 3,  useDingbats=FALSE )
print(s_fig) 
dev.off()


