############################################
### plotting methods for protrider
############################################


# pacman::p_load(data.table, magrittr, matrixStats, ggplot2, parallel, ggpubr,plotly, 
#                keras, ROCR, pheatmap, pcaMethods, patchwork, SummarizedExperiment,
#                VennDiagram, gridExtra, MASS, robustbase, RColorBrewer, plyr)

pacman::p_load(ggplot2, plotly, pheatmap, VennDiagram, gridExtra, RColorBrewer, cowplot, gplots)



### small preview of matrix
pre = function(df) {
  rows = min(nrow(df),6)
  cols = min(ncol(df),6)
  df[1:rows,1:cols]
}

### sorts second matrix row and cols like the first
sort_plot_input = function(x1,x2) {
    if(!is.matrix(x2)){
        return(x2)
    }
    x2_row_sort = match(rownames(x1), rownames(x2))
    row_x2 = x2[x2_row_sort,]  
    
    x2_col_sort = match(colnames(x1), colnames(row_x2))
    col_row_sorted = x2[,x2_col_sort]
    
    return(col_row_sorted)
}



############################################
#' # QQ-plots

qq_sampleProt = function(prot, pval, outlier_pos="", main="p-values QQ-plot", conf.alpha=0.05){
    # if(outlier_pos =="") {
    #     title_text = ''
    # } else {
    #     title_text = paste0('\n[marked: ',outlier_pos,']')
    # }
    
    obsP = pval[prot,]
    obsP = sort(na.omit(obsP))
    expP = ppoints(obsP)
    
    df = data.frame(obsP = -log10(obsP), expP = -log10(expP), is_marked = names(obsP) == outlier_pos, rank=1:length(obsP))
    df = df[order(df$expP),]
    p = ggplot(df, aes(x = expP, y=obsP, col=is_marked)) +
        geom_point(color=c("black", "firebrick")[df$is_marked + 1]) +
        theme_bw() + theme(legend.position="none", plot.title = element_text(hjust = 0.5)) +
        ggtitle(main) + #xlab("expected p-value [-log10]") + ylab("observed p-value [-log10]")
        xlab(expression(-log[10]~"(expected p-value)")) +
        ylab(expression(-log[10]~"(observed p-value)"))
        
    if(is.numeric(conf.alpha)) {
        df$upper = -log10(qbeta(conf.alpha/2, df$rank, max(df$rank) - df$rank))
        df$lower = -log10(qbeta(1-conf.alpha/2, df$rank, max(df$rank) - df$rank))
        p = p + geom_ribbon(data=df, aes(x=expP, ymin=lower, ymax=upper, text=NULL), alpha=0.2, color="gray")
        }
    p = p + geom_abline(intercept=0, slope=1, col="orange")
    p
  }


### fast version of qqplots
qq_globalProt = function(pval, outlier_pos = NULL, main="global QQ plot"){
    obsP = c(unlist(-log10(pval)))
    if(!is.null(outlier_pos)){
        outlier_pos[is.na(pval)] = NA
        outP = c(unlist(outlier_pos))
        ggplot_dt = data.table(obsP = obsP, inj_outlier=outP)
    } else{
        ggplot_dt = data.table(obsP = obsP)
    }
    ggplot_dt = na.omit(ggplot_dt)
    ggplot_dt = ggplot_dt[order(ggplot_dt$obsP, decreasing = T),]
    ggplot_dt$expP = c(unlist(-log10(ppoints(ggplot_dt$obsP))))
    ggplot_dt[ggplot_dt>20] = 20

    pl = ggplot(ggplot_dt, aes(x=expP,y=obsP)) +
        stat_binhex(aes(fill=log10(..count..))) +
        geom_abline(intercept = 0, slope = 1, color='orange') +
        ggtitle(main) + ylim(c(0,20)) + ylab(expression(-log[10]~"(observed p-value)")) + xlab(expression(-log[10]~"(expected p-value)")) +
        theme(legend.position = "none")

    if(!is.null(outlier_pos)){
        pl + geom_point(data = subset(ggplot_dt,inj_outlier !=0), mapping = aes(x = expP, y = obsP), color="red")
    }
    else{
        pl
    }
}

### slower version of qq-plot but more detailed
qqplot_global_detailed = function(pval, main="",  conf.alpha=0.05){
  obsP = as.vector(pval)
  obsP = sort(na.omit(obsP))
  expP = ppoints(obsP)
  
  df = data.frame(obsP = -log10(obsP), expP = -log10(expP), rank=1:length(obsP))
  df = df[order(df$expP),]
  p = ggplot(df, aes(x = expP, y=obsP)) +
    geom_point() +
    theme_bw() + theme(text = element_text(size=15)) +  # ggtitle(main) +
    xlab(expression(-log[10]~"(expected p-value)")) + ylab(expression(-log[10]~"(observed p-value)"))
  
  #theme(legend.position="none", text = element_text(size=15) )
  if(is.numeric(conf.alpha)) {
    df$upper = -log10(qbeta(conf.alpha/2, df$rank, max(df$rank) - df$rank))
    df$lower = -log10(qbeta(1-conf.alpha/2, df$rank, max(df$rank) - df$rank))
    p = p + geom_ribbon(data=df, aes(x=expP, ymin=lower, ymax=upper, text=NULL), alpha=0.2, color="gray")
  }
  p = p + geom_abline(intercept=0, slope=1, col="orange")
  p
}






############################################
#' # helper methods for plotting
plot_binhex = function(x, y, title="", axis_limit=20, center=TRUE, xylabel=c("plot_X", "plot_Y")){
  y = sort_plot_input(x,y)
  axis_min = -axis_limit
  # if(axis_limit>30) {
  #     axis_min = max(axis_min, 0)
  # }
  ggplot_dt = data.table(plot_X = c(unlist(x)), plot_Y = c(unlist(y)) )
  pl = ggplot(ggplot_dt, aes(x=plot_X,y=plot_Y)) +
    stat_binhex(aes(fill=log10(..count..))) +
    geom_abline(intercept = 0, slope = 1, color='orange') + ggtitle(title) +
      xlab(xylabel[1]) + ylab(xylabel[2]) + theme(legend.position = "none")
  
  if(center) {
    pl + xlim(axis_min, axis_limit) + ylim(axis_min, axis_limit) 
  } else {
    pl
  }
}


plot_heatscatter = function(se, x_ax, y_ax, main, xylabel=c("plot_X", "plot_Y"), axis_limit=7) {
  # x_val = 
  df = data.frame(plot_X = as.vector(assay(se, x_ax)), plot_Y = as.vector(assay(se, y_ax)), is_na = as.vector(assays(se)$X_na ))
  df = subset(df, is_na!=TRUE)
  LSD::heatscatter(df$plot_X, df$plot_Y,  main=main,
                   xlim=c(-axis_limit,axis_limit), ylim=c(-axis_limit,axis_limit),
                   xlab=xylabel[1], ylab=xylabel[2])
  abline(0,1, col="orange")
}




plot_pvalues_binhex = function(pValues_run1, pValues_run2, outliers_run1=NULL, title="") {
    pValues_run2 = sort_plot_input(pValues_run1, pValues_run2)
    
    ggplot_df = data.frame(run1=-log10(as.vector(pValues_run1)), run2=-log10(as.vector(pValues_run2)))
    if(!is.null(outliers_run1)){
        ggplot_df$outliers=as.vector(outliers_run1)
    } 
    ggplot_df = na.omit(ggplot_df)
    ggplot_df = ggplot_df[!is.infinite(rowSums(ggplot_df)),]  # inf values if pvalues too small
    
    run_cor = cor(c(ggplot_df$run1), c(ggplot_df$run2) ,use = "complete.obs")  # nother one is due to there being one of the values being constant
    ggplot_df[ggplot_df>25] = 25
    #out_plot = ggplot(ggplot_df, aes(y = run1, x=run2, color=density)) + scale_color_viridis() + + #scale_x_log10() + scale_y_log10()+ 
    pl = ggplot(ggplot_df, aes(y = run1, x=run2)) +
        stat_binhex(aes(fill=log10(..count..))) + ggtitle(paste0(title, " [cor:",round(run_cor,3),"]")) + xlab("run2") + ylab("run1") + ylim(c(0,25))+xlim(c(0,25))+
        geom_hline(yintercept=5, color = "red") + geom_vline(xintercept=5, color = "darkred") + geom_abline(intercept = 0, slope = 1, color='orange')
    if(!is.null(outliers_run1)){
        pl + geom_point(data = subset(ggplot_df, outliers !=0), mapping = aes(y = run1, x=run2), color="red")
    }
    else {
        pl
    }
}

plot_pvalues_scatter = function(pValues_run1, pValues_run2, outlier_pos, title="") {
    pvalue_df = data.frame(run1=-log10(as.vector(pValues_run1)), run2=-log10(as.vector(pValues_run2)), is_outlier=outlier_pos)
    pvalue_df = na.omit(pvalue_df)
    # run_cor = cor(pvalue_df$run1, pvalue_df$run2 ,use = "complete.obs")
    pvalue_df[pvalue_df>25] = 25
    
    LSD::heatscatter(pvalue_df$run1, pvalue_df$run2,  main=title, 
                     xlab = "run 1 p-value [-log10]", ylab="run 2 p-value [-log10]",
                     xlim=c(0,26), ylim=c(0,26), cor=TRUE, 
                     colpal="standard"       )
    abline(0,1, col="orange")
}


### requires: data_handling/get_prot_sample_list() as innput
plot_venn = function(df1, df2, sig_level=0.1, fold_change=1.5, df_label=c("df1","df2"), main="sample_protein overlap outiers") {
    df1_samples = unlist(as.vector(subset(df1, adj_pvalue<sig_level & abs(log2fc) > log2(fold_change) )[,'sample_prot']) )
    df2_samples = unlist(as.vector(subset(df2, adj_pvalue<sig_level & abs(log2fc) > log2(fold_change) )[,'sample_prot']))
    
    # df1_samples = as.vector(subset(df1,is_outlier==TRUE) [,'sample_prot'])
    # df2_samples = as.vector(subset(df2,is_outlier==TRUE) [,'sample_prot'])
    
    
    df_intersect = length(intersect(df1_samples,df2_samples))
    grid.newpage()
    g = draw.pairwise.venn(area1 = length(df1_samples), area2 = length(df2_samples), cross.area = df_intersect, category = df_label,  fill = c("green", "gray"), 
                            alpha = rep(0.4, 2), cat.pos = c(337, 25),
                            lty = rep("blank", 2), cex = 1.5, cat.cex = 1.5, fontfamily=rep("Arial",3), cat.fontfamily = rep("Arial", 2)
    )#euler.d = F, sep.dist = 0.03, rotation.degree = 45)
    if(is.null(main)) {
        g
    } else {
        # grid.arrange(gTree(children=g), top=paste0(main)) #,' [adj_pval<',sig_level,' & fc>',fold_change,']'))
        grid.arrange(gTree(children=g), top=textGrob(main, gp=gpar(fontsize=15,font=8))) #,' [adj_pval<',sig_level,' & fc>',fold_change,']'))
    }
}




### requires: data_handling/get_prot_sample_list() with outliers as innput
plot_venn_outlier = function(df1, df2, sig_level=0.1, fold_change=1.5) {
    df1_samples = as.vector(subset(df1, adj_pvalue<sig_level & abs(log2fc) > log2(fold_change) )[,'sample_prot'])
    df2_samples = as.vector(subset(df2, adj_pvalue<sig_level & abs(log2fc) > log2(fold_change) )[,'sample_prot'])
    
    # df1_samples = as.vector(subset(df1, is_outlier==TRUE )[,'sample_prot'])
    # df2_samples = as.vector(subset(df2, is_outlier==TRUE )[,'sample_prot'])
    
    df3_samples = as.vector(subset(df1, outliers!=0)[,'sample_prot'])
    
    in12 = length(intersect(df1_samples,df2_samples))
    in23 = length(intersect(df2_samples,df3_samples))
    in13 = length(intersect(df1_samples,df3_samples))
    in123 = length(intersect(intersect(df1_samples,df2_samples),df3_samples))
    
    grid.newpage()
    g = draw.triple.venn(
        area1 = length(df1_samples), area2 = length(df2_samples), area3 = length(df3_samples),
        n12 = in12,
        n23 = in23,
        n13 = in13,
        n123 = in123,
        category = c("run1", "run2", "injected_outliers"),
        fill = c("light blue", "green", "red"),
        lty = "blank", cex = 1.5, cat.cex = 1.5
    )
    grid.arrange(gTree(children=g), top=paste0('sample_protein overlap   [adj_pval<',sig_level,' & fc>',fold_change,']'))

}



### TODO CHANGE <0.1 to outlier_called of se
### plots single protein data for prokisch format
plot_single_prot = function(se, protein, marked_sample=NULL){
    par(mfrow=c(1,3))
    # plot(assays(se)$X_raw[protein,], log="y", col=1+(assays(se)$X_out_called), main=paste0("raw data ",protein), xlab="samples")
    # if(!is.null(marked_sample)) {   points(which(names(assays(se)$X_raw[protein,]) %in% marked_sample), assays(se)$X_raw[protein,marked_sample], col="green", pch=16)  }
    
    plot(assays(se)$X[protein,], col=1+(assays(se)$X_out_called), main=paste0(protein,"  centered data"), xlab="samples")
    if(!is.null(marked_sample)) {   points(which(names(assays(se)$X[protein,]) %in% marked_sample), assays(se)$X[protein,marked_sample], col="green", pch=16)  }
    
    plot(assays(se)$X_pred[protein,], col=1+(assays(se)$X_out_called), main="prediction", xlab="samples")
    if(!is.null(marked_sample)) {   points(which(names(assays(se)$X_pred[protein,]) %in% marked_sample), assays(se)$X_pred[protein,marked_sample], col="green", pch=16)  }
    
    
    plot(assays(se)$X_pred[protein,],assays(se)$X[protein,], col=1+(assays(se)$X_out_called), ylab="centered", xlab="prediction", main="red: significant, blue: marked_sample")
    if(!is.null(marked_sample)) {   points(assays(se)$X_pred[protein,marked_sample], assays(se)$X[protein,marked_sample], col="green", pch=16)  }
    abline(0,1, col="blue")
    par(mfrow=c(1,1))
}







# sig_level = 0.1
# df1 = se_pval_list
# df2 = se2_pval_list
# df_label=c("df1","df2")

plot_foldchange_boxplot3 = function(df1, df2, sig_level=0.1,fold_change=1.5, df_label=c("df1","df2") ) {
    df1_samples = as.vector(subset(df1, adj_pvalue<sig_level & abs(log2fc) > log2(fold_change) )[,'sample_prot'])
    df2_samples = as.vector(subset(df2, adj_pvalue<sig_level & abs(log2fc) > log2(fold_change) )[,'sample_prot'])
    df_intersect =intersect(df1_samples,df2_samples)
    
    df1_samples = df1_samples[!df1_samples %in% df_intersect]
    df2_samples = df2_samples[!df2_samples %in% df_intersect]
    
    df1_box = subset(df1, sample_prot %in% df1_samples)
    if(nrow(df1_box)>0) {
        df1_box$group = df_label[1]
    }
    df2_box = subset(df2, sample_prot %in% df2_samples)
    if(nrow(df2_box)>0) {
        df2_box$group = df_label[2]
    }
    df_intersect_box = subset(rbind(df1,df2), sample_prot %in% df_intersect)
    df_intersect_box$group = "intersect"
    
    plot_df = rbind(df1_box, df2_box, df_intersect_box)
    
    # Perform pairwise comparisons
    compare_means(log2fc~group,data=plot_df, method = "wilcox.test", paired = FALSE)
    
    # add sig comparisons 
    my_comparisons = list( c(df_label[1], df_label[2]), c(df_label[2], "intersect"), c(df_label[1], "intersect") )
    plot_df$group = as.factor(plot_df$group)
    
    g = ggplot(plot_df, aes(x =group, y = log2fc)) + geom_violin() +
        stat_compare_means(comparisons = my_comparisons) +
        ggtitle("foldchange of venn regions [wilcox unpaired]") +
        geom_hline(aes(yintercept = 0), linetype = 2)
    
    return(list(plot_df = plot_df, plot = g))
}




# protein="P41567"
# marked_sample ="P102386"
    

### plots single protein for 2 runs, marked sample from se
plot_single_prot_injOutliers = function(se, se2, protein, marked_sample=NULL){
    
    box_df = as.matrix(DataFrame("se"=assays(se)$X_pred[protein,] - assays(se)$X[protein,], "se2"=assays(se2)$X_pred[protein,] - assays(se2)$X[protein,], "pred_se_minus_se2" = assays(se)$X_pred[protein,] -assays(se2)$X_pred[protein,],"se_adj_pval" = assays(se)$X_pval_adj[protein,]))
    pred_xlim = c(min(assays(se)$X_pred[protein,],assays(se2)$X_pred[protein,], na.rm = T) -0.2, max(assays(se)$X_pred[protein,],assays(se2)$X_pred[protein,], na.rm = T)+0.2 )
    pred_ylim = c(min(assays(se)$X[protein,],assays(se2)$X[protein,],na.rm = T)-0.2, max(assays(se)$X[protein,], assays(se2)$X[protein,],na.rm = T)+ 0.2 )
    
    cumsum(colCounts(as.matrix(colData(se))==1))
    
    par(mfrow=c(2,3))
    plot(assays(se)$X_raw[protein,], log="y", col=1+(assays(se)$X_pval_adj[protein,]<0.1), main=paste0("raw data ",protein), xlab="samples", pch=16)
    if(!is.null(marked_sample)) {   points(which(names(assays(se)$X_raw[protein,]) %in% marked_sample), assays(se)$X_raw[protein,marked_sample], col="green", pch=0)  }
    points(which(assays(se)$X_out_pos[protein,] != 0), assays(se)$X_raw[protein ,assays(se)$X_out_pos[protein,] != 0 ], col="blue", pch=4)
    abline(v = cumsum(colCounts(as.matrix(colData(se))==1)), col="grey")
    
    
    plot(assays(se)$X[protein,], col=1+(assays(se)$X_pval_adj[protein,]<0.1), xlab="samples", main="red: significant, blue:injected_outlier, green: marked_sample", pch=16)
    if(!is.null(marked_sample)) {   points(which(names(assays(se)$X[protein,]) %in% marked_sample), assays(se)$X[protein,marked_sample], col="green", pch=0)  }
    points(which(assays(se)$X_out_pos[protein,] != 0), assays(se)$X[protein ,assays(se)$X_out_pos[protein,] != 0 ], col="blue", pch=4)
    abline(0,0, col="grey")
    abline(v = cumsum(colCounts(as.matrix(colData(se))==1)), col="grey")
    

    plot(assays(se)$X_pred[protein,], col=1+(assays(se)$X_pval_adj[protein,]<0.1), xlab="samples", main="red: significant, blue:injected_outlier, green: marked_sample", pch=16)
    if(!is.null(marked_sample)) {   points(which(names(assays(se)$X[protein,]) %in% marked_sample), assays(se)$X[protein,marked_sample], col="green", pch=0)  }
    points(which(assays(se)$X_out_pos[protein,] != 0), assays(se)$X[protein ,assays(se)$X_out_pos[protein,] != 0 ], col="blue", pch=4)
    abline(0,0, col="grey")
    abline(v = cumsum(colCounts(as.matrix(colData(se))==1)), col="grey")

#    boxplot(box_df[,1:3], ylim=pred_ylim, main="foldchange to centered")
#    abline(0,0, col="grey")
    
    plot(assays(se)$X_pred[protein,], assays(se)$X[protein,], col=1+(assays(se)$X_pval_adj[protein,]<0.1), ylab="centered", xlab="prediction", main="se pred", ylim=pred_ylim, xlim=pred_xlim)
    if(!is.null(marked_sample)) {   points(assays(se)$X_pred[protein,marked_sample], assays(se)$X[protein,marked_sample], col="green", pch=0)  }
    points(assays(se)$X_pred[protein ,assays(se)$X_out_pos[protein,] != 0 ], assays(se)$X[protein ,assays(se)$X_out_pos[protein,] != 0 ], col="blue", pch=4)
    abline(0,1, col="grey")
    
    plot(assays(se2)$X_pred[protein,], assays(se2)$X[protein,], col=1+(assays(se2)$X_pval_adj[protein,]<0.1), ylab="centered", xlab="prediction", main="se2 pred", ylim=pred_ylim, xlim=pred_xlim)
    if(!is.null(marked_sample)) {   points(assays(se2)$X_pred[protein,marked_sample], assays(se2)$X[protein,marked_sample], col="green", pch=0)  }
    points(assays(se2)$X_pred[protein ,assays(se2)$X_out_pos[protein,] != 0 ], assays(se2)$X[protein ,assays(se2)$X_out_pos[protein,] != 0 ], col="blue", pch=4)
    abline(0,1, col="grey")
    
    
    plot(box_df[,"pred_se_minus_se2"], main="prediction difference" ) #, col=1+(assays(se)$X_pval_adj[protein,]<0.1)) #col=1+(box_df[,"se_adj_pval"]<0.1) )
    if(!is.null(marked_sample)) {   points(which(names(assays(se)$X[protein,]) %in% marked_sample), box_df[,"pred_se_minus_se2"][which(names(assays(se)$X[protein,]) %in% marked_sample)], col="green", pch=0)  }
    #points(which(box_df[,"se_adj_pval"]<0.1), box_df[which(box_df[,"se_adj_pval"]<0.1),"pred_se_minus_se2"], col="red")
    #points(which(assays(se)$X_pval_adj[protein,]<0.1), box_df[which(assays(se)$X_pval_adj[protein,]<0.1),"pred_se_minus_se2"], col="red")
    points(which(assays(se)$X_out_pos[protein,] != 0), box_df[,"pred_se_minus_se2"][which(assays(se)$X_out_pos[protein,] != 0)], col="blue", pch=4)
    
    abline(0,0, col="grey")
    
    par(mfrow=c(1,1))
}


plot_heatmap = function(se, plot_matrix, main="", confounders_used=c("PROTEOMICS_BATCH","gender","INSTRUMENT"), sample_id_col = "proteome_ID", ... ){
    if(is.null(metadata(se)$batches)) {
        anno_col = get_batches_melted(colData(se))
    } else {
        
        anno_col = subset(metadata(se)$batches, USE_FOR_PROTEOMICS_PAPER == TRUE | USE_FOR_PROTEOMICS_PAPER=="True")
        anno_col = anno_col[, c(sample_id_col, confounders_used)]
        rownames(anno_col) = as.character(anno_col[[sample_id_col]])
        anno_col[sample_id_col] = NULL
        
        ### transform to factors
        for(con in confounders_used) {
          anno_col[con] = as.factor(anno_col[[con]])
        }

      }
  
    plot_matrix[assay(se, "X_na")] = NA  
  
    p = pheatmap(cor(plot_matrix, use="complete.obs"), main=main, 
             breaks = seq(-1,1, length.out = 100),na_col = "grey", annotation_col = anno_col, 
             color=colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100), ...)
    p
}








plot_aberrant_sample = function(se, main, padjCutoff=0.05, zScoreCutoff=0,
                                  outlierRatio=0.001,
                                  # col=brewer.pal(3, 'Dark2')[c(1,2)],
                                 col=c('orange','blue'),
                                    yadjust=c(1.2, 1.2), 
                                  labLine=c(3.5, 3), ylab="# aberrantly expressed proteins", 
                                  labCex=par()$cex, ymax=NULL, ...){
    
    if(missing(main)){
        main = 'aberrant proteins per sample'
    }
    
    if( class(se) =="SummarizedExperiment") {
        count_vector = sort(get_aberrant(se, by="sample", padjCutoff=padjCutoff) )
        se_length = length(se)
    } else {
        count_vector = sort(table(subset(se, PROT_PADJ<padjCutoff)$SAMPLEID))  ## only for limma_raw input !!
        # count_vector = sort(table(subset(se, signif==TRUE)$SAMPLEID))  ## only for limma_raw input !!
        se_length = length(unique(se$UNIPROT_ID))
        }
    
                         
    ylim = c(0.4, max(1, count_vector)*1.1)
    if(!is.null(ymax)){
        ylim[2] = ymax
    }
    replace_zero_unknown = 0.5
    ticks= c(replace_zero_unknown, signif(10^seq(
        from=0, to=round(log10(max(1, count_vector))), by=1/3), 1))
    
    labels_for_ticks = sub(replace_zero_unknown, '0', as.character(ticks))
    
    bp= barplot2(
        replace(count_vector, count_vector==0, replace_zero_unknown),
        log='y', ylim=ylim, names.arg='', xlab='', plot.grid=TRUE, 
        grid.col='lightgray', ylab='', yaxt='n', border=NA, xpd=TRUE,
        col=col[(!count_vector <= max(1, se_length*outlierRatio)) + 1],
        main=main)
    
    n_names = floor(length(count_vector)/20)
    xnames= seq_len(n_names*20)
    axis(side=1, at= c(0,bp[xnames,]), labels= c(0,xnames))
    axis(side=2, at=ticks, labels= labels_for_ticks, ylog=TRUE, las=2)
    
    # labels
    mtext('sample rank', side=1, line=labLine[1], cex=labCex)
    mtext(ylab, side=2, line=labLine[2], cex=labCex)
    
    # legend and lines
    hlines = c(median=ifelse(median(count_vector)==0, replace_zero_unknown,
                             median(count_vector)) , Quantile90=quantile(
                                 count_vector,0.9, names=FALSE))
    color_hline= c('black','black')
    abline(h= hlines, col=color_hline)
    text(x=c(1,1), y= hlines*yadjust, col=color_hline, adj=0,
         labels=c('median', expression(90^th ~ 'percentile')))
    
    box()
}




### proteins on y and samples on x axis
plot_protein_scatter = function(se, plot_table, sample_id, protein_id, main=NULL, legend=TRUE) {
    index = which(rownames(se) == protein_id)
    # plot_table = assays(se)$X
    plot_table[assays(se)$X_na] = NA
    batches = get_batches_melted(colData(se))
    plot_obj = list("prot_data"=as.data.frame(t(plot_table)), "batches"=batches)
    
    if(is.null(main)) {
        title_text = paste0(protein_id," [marked sample: ",sample_id,"]")
    } else {
        title_text = main
    }
    
    ggplot_df = data.frame(index= 1:nrow(plot_obj$prot_data), intensity=plot_obj$prot_data[[protein_id]], sample =rownames(plot_obj$prot_data), batch =plot_obj$batches$batches)
    p = ggplot(ggplot_df, aes(x=index, y=intensity, colour = batch, name=sample)) + geom_point() + ggtitle(title_text) +
        #scale_color_brewer(palette="Set1")
        scale_colour_manual(values=c("#89C5DA", "#DA5724", "#74D944", "#CE50CA", "#3F4921", "#C0717C", "#CBD588", "#5F7FC7",
                                     "#673770", "#D3D93E", "#38333E", "#508578", "#D7C1B1", "#689030", "#AD6F3B", "#CD9BCD",
                                     "#D14285", "#6DDE88", "#652926", "#7FDCC0", "#C84248", "#8569D5", "#5E738F", "#D1A33D",
                                     "#8A7C64", "#599861")) +
        theme_bw() + 
        geom_point(data = subset(ggplot_df, sample == sample_id), mapping = aes(x=index, y=intensity), color="orange", size=4, alpha=0.3)  +
        xlab("samples") + ylab(expression(log[2]~"(intensity)")) #ylab("intensity [log2]")
    
    if(legend) {
        p = p + theme(legend.position="bottom")
    } else {
        p = p + theme(legend.position="None")
    }
    
    return(p)
}

plot_protein_scatter_all = function(se, sample_id, protein_id) {
    p_raw = plot_protein_scatter(se, assays(se)$X + assays(se)$X_mean, sample_id=sample_id, protein_id=protein_id, main="observed", legend=FALSE)
    p_pred = plot_protein_scatter(se, assays(se)$X_pred + assays(se)$X_mean, sample_id=sample_id, protein_id=protein_id, main="expected", legend=FALSE)
    p_cont = plot_protein_scatter(se, assays(se)$X - assays(se)$X_pred, sample_id=sample_id, protein_id=protein_id, main="contrast", legend=FALSE)
    p = plot_grid(p_raw, p_pred, p_cont, labels = "AUTO", ncol=3)
    
    title = ggdraw() + draw_label(paste0('protein ',protein_id,' with marked sample ',sample_id), fontface='bold')
    plot_grid(title, p, ncol=1, rel_heights=c(0.1, 1))

}



### prot_list are intensity row of a protein
plot_protein_ranked = function(prot_row, sample_id="", main="", ylab=expression(log[2]~"(protein intensity)")) {
    prot_row_rm = prot_row[prot_row!=0]
    prot_row_rm = sort(prot_row_rm[is.finite(prot_row_rm)])
    df = data.frame(prot_intensity = prot_row_rm, is_marked = names(prot_row_rm)==sample_id, sample_rank = 1:length(prot_row_rm))
    p = ggplot(df, aes(x=sample_rank, y=prot_intensity, col=is_marked))+
        geom_point(alpha=ifelse(df$is_marked, 1, 0.5),
                   color=c("gray70", "firebrick")[df$is_marked + 1]) +
        theme_bw() +
        theme(legend.position="none", plot.title = element_text(hjust = 0.5)) +
        ggtitle(main) + xlab('sample rank') + ylab(ylab)
    p
}

# plot_protein_ranked((assays(se)$X + assays(se)$X_mean)[1,], sample_id = "P103216")


### ranked raw - protrider - limma intensity
plot_protein_ranked_comparison = function(se, limma_table_raw, protein_id, sample_id, main=""){
    library(cowplot)
    
    ### data handling
    prot_before = (assays(se)$ X+ assays(se)$X_mean)[protein_id,]
    
    limma_prot = subset(limma_table_raw, UNIPROT_ID==protein_id)
    prot_limma = limma_prot$NORM_LOG2_LFQ
    names(prot_limma) = limma_prot$SAMPLEID
    
    ### plot
    raw = plot_protein_ranked(prot_before, sample_id, main="raw intensities")
    protrider_norm = plot_protein_ranked( prot_before - assays(se)$X_pred[protein_id,], sample_id, main="PROTRIDER normalized")
    limma_norm = plot_protein_ranked( prot_limma, sample_id, main="LIMMA normalized")
    p = plot_grid(raw, protrider_norm, limma_norm, labels = "AUTO", ncol=3)

    title = ggdraw() + draw_label(paste0('protein ',protein_id,' with marked sample ',sample_id), fontface='bold')
    plot_grid(title, p, ncol=1, rel_heights=c(0.1, 1))
}

plot_protein_ranked_comparison_short = function(se, protein_id, sample_id, main=""){
    prot_before = (assays(se)$ X+ assays(se)$X_mean)[protein_id,]
    
    raw = plot_protein_ranked(prot_before, sample_id, main="raw intensities")
    protrider_norm = plot_protein_ranked( prot_before - assays(se)$X_pred[protein_id,], sample_id, main="PROTRIDER normalized")
    p = plot_grid(raw, protrider_norm, labels = "", ncol=2)
    p
}

# plot_protein_ranked_comparison(se, limma_raw, protein_id = "A0AVF1", sample_id = "P103216")


### adapted from outrider, ods is se object
plotVolcano = function(ods, sampleID, main, x_value=c("fc","z"), protein_marked="", padjCutoff=0.05, zScoreCutoff=0,
                        pch=16, basePlot=TRUE, col=c("gray70", "firebrick", "blue")) {
    
    dt = data.table(
        GENE_ID   = rownames(ods),
        pValue    = assays(ods)$X_pval[,sampleID],
        padjust   = assays(ods)$X_pval_adj[,sampleID],
        zScore    = assays(ods)$X_z[,sampleID],
        fc    = 2^assays(ods)$X_log2fc[,sampleID],
        normCts   = assays(ods)$X[,sampleID],
        medianCts = rowMedians(assays(ods)$X),
        expRank   = apply( assays(ods)$X, 2, rank)[,sampleID] ,
        aberrant  = assays(ods)$X_out_called[,sampleID] ,
        color=col[1])
    dt[aberrant == TRUE, color:=col[2]]
    dt[aberrant == TRUE, color:=col[2]]
    
    dt[GENE_ID == protein_marked, color:=col[3]]
    
    x_plot = if (x_value=="fc") "fc" else "zScore" 
    
    # remove the NAs from the zScores for plotting
    dt[is.na(zScore),zScore:=0]
    # dt[is.na(log2fc),log2fc:=0]
    
    if(isTRUE(basePlot)){
        p = ggplot(dt, aes(x=get(x_plot), y=-log10(pValue)) )+
            geom_point(color=dt$color )+
            theme_bw() +
            theme(legend.position="none", plot.title = element_text(hjust = 0.5)) +
            ggtitle(main) + xlab(x_value) + ylab(expression(-log[10]~"(p-value)"))
        return(p)
        
    }
    plot_ly(
        data=dt,
        x=~get(x_plot),
        y=~-log10(pValue),
        #y =~pValue,
        type="scatter",
        mode="markers",
        marker = list(
            color=~color
        ),
        text=~paste0(
            "Gene ID: ", GENE_ID,
            "<br>Sample ID: ", sampleID,
            "<br>Median normcount: ", round(medianCts, 2),
            "<br>normcount: ", round(normCts, 2),
            "<br>expression rank: ", as.integer(expRank),
            "<br>nominal P-value: ", signif(pValue,3),
            "<br>adj. P-value: ", signif(padjust,3),
            "<br>Z-score: ", signif(zScore,2),
            "<br>fc: ", signif(fc,2)
        )
    ) %>% layout(title=main, 
                 yaxis=list(title="-log<sub>10</sub>(<i>P</i>-value)"))
}
# plotVolcano(se, sample_id, main="", protein_marked=protein_id, col=c("gray70", "gold2", "firebrick"))




plot_observed_expected = function(se, sample_id, protein_id, main="observed against expected intensities") {
    X = (assays(se)$X + assays(se)$X_mean)
    X[assays(se)$X_na] = NA
    X_pred = (assays(se)$X_pred + assays(se)$X_mean)
    X_pred[assays(se)$X_na] = NA
    
    plot_df = data.frame(sample_id=colnames(se), X = X[protein_id,], X_pred = X_pred[protein_id,])
    plot_df$is_marked = plot_df$sample_id == sample_id
    p = ggplot( plot_df, aes(x=X_pred, y=X, col=is_marked)) +
        geom_point(color=c("black", "firebrick")[plot_df$is_marked + 1]) +
        theme_bw() + theme(legend.position="none", plot.title = element_text(hjust = 0.5)) +
        ggtitle(main ) +
        # xlab('expected intensity [log2]') +
        # ylab('observed intensity [log2]') +
        xlab(expression(log[2]~"(expected intensity)")) + 
        ylab(expression(log[2]~"(observed intensity)")) +
        geom_abline(intercept = 0, slope = 1, color='orange')
    p
}
    
    

### scatter proteins with control and outliers confirmed / found marked
plot_control_outlier_check_scatter = function(se_control, prot_id, sa_prot, title_text = "", sample_id_col = "proteome_ID") {
  index = which(rownames(se_control) == prot_id)
  plot_table = assays(se_control)$X
  plot_table[assays(se_control)$X_na] = NA
  batches = get_batches_melted(colData(se_control))
  plot_obj = list("prot_data"=as.data.frame(t(plot_table)), "batches"=batches)
  
  ggplot_df = data.frame(index= 1:nrow(plot_obj$prot_data), intensity=plot_obj$prot_data[[prot_id]],
                         sample =rownames(plot_obj$prot_data), batch =plot_obj$batches$batches)
  
  outlier_df = data.frame(sample = colnames(se_control), is_outlier = assays(se_control)$X_out_called[prot_id,])
  ggplot_df = merge(ggplot_df, outlier_df)
  ggplot_df$confirmed = sapply(ggplot_df$sample, function(x) { x %in% subset(sa_prot, protein_id == prot_id)[sample_id_col] })
  ggplot_df$is_control = grepl("PS", ggplot_df$sample)
  
  p = ggplot(ggplot_df, aes(x=index, y=intensity, colour = batch, name=sample))  + ggtitle(title_text) +
    #scale_color_brewer(palette="Set1")
    scale_colour_manual(values=c("#89C5DA", "#DA5724", "#74D944", "#CE50CA", "#3F4921", "#C0717C", "#CBD588", "#5F7FC7",
                                 "#673770", "#D3D93E", "#38333E", "#508578", "#D7C1B1", "#689030", "#AD6F3B", "#CD9BCD",
                                 "#D14285", "#6DDE88", "#652926", "#7FDCC0", "#C84248", "#8569D5", "#5E738F", "#D1A33D",
                                 "#8A7C64", "#599861")) +
    theme_bw() + 
    geom_point(data = subset(ggplot_df, confirmed==TRUE), mapping = aes(x=index, y=intensity), color="red", size=5, pch=1)  +
    geom_point(data = subset(ggplot_df, is_outlier==TRUE), mapping = aes(x=index, y=intensity), color="orange", size=4, alpha=0.5)  +
    geom_point(data = subset(ggplot_df, is_control==TRUE), pch=1 )  +
    geom_point(data = subset(ggplot_df, is_control==FALSE), pch=19)  +
    xlab("samples") + theme(legend.position="None") + ylab(expression(log[2]~"(intensity)"))
  p
}

    

### diverse mix of multiple plots
plot_full_check = function(se, sample_id, protein_id){
  p1 = plot_protein_ranked_comparison_short(se, sample_id=sample_id, protein_id=protein_id)
  p2 = qq_sampleProt(prot=protein_id, pval=assays(se)$X_pval, outlier_pos=sample_id, main="QQ-plot p-values")
  p3 = plot_protein_ranked(log10(assays(se)$X_pval[protein_id,]), sample_id = sample_id, main="ranked p-values", ylab=expression(-log[10]~"(p-value)") )
  p4 = plot_protein_scatter(se, assays(se)$X + assays(se)$X_mean, sample_id=sample_id, protein_id=protein_id, main=paste0("protein ",protein_id), legend=FALSE)
  p5 = plotVolcano(se, sample_id, main=sample_id, protein_marked=protein_id,x_value="z", col=c("gray70", "gold2", "firebrick"))
  # p5 = plotVolcano(se, sample_id, main=sample_id, protein_marked=protein_id,x_value="fc", col=c("gray70", "gold2", "firebrick"))
  p6 = plot_observed_expected(se, sample_id, protein_id, main=protein_id)
  plot_grid(p1,p2,p3,p4,p5,p6, labels = "AUTO", nrow=2, rel_widths=c(4,2,2), scale = 0.9)
}



