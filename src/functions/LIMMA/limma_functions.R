suppressPackageStartupMessages({
    stringsAsFactors=F
    library(readxl)
    library(readr)
    library(openxlsx)
    require(data.table)
    library(dplyr)
    library(tidyr)
    
    library(DESeq2)
    library(pracma)
    library(scales)
    require(magrittr)
    require(limma)
    
    library(biomaRt)
    library(org.Hs.eg.db)
    library(ggrepel)
    
    library(ggplot2)
    library(gplots)
    library(ggthemes)
    library(ggpubr)
    require(pheatmap)
    require(RColorBrewer)
})





## remove zeroblocks
rm_zeroblocks <- function(data, blocks) {
    rm_rows <- NULL
    for(i in 0:(blocks-1)) {
        tmp_rm_rows <- rowSums(!data[, (1+i*10):(10+i*10)])==10
        if (i == 0) {
            rm_rows <- tmp_rm_rows
        } else {
            rm_rows <- rm_rows | tmp_rm_rows
        }
    }
    data.nozeroblocks <- data[!rm_rows, ]
    return(data.nozeroblocks)
}


#histogram function
histnorm3 <- function(COUNT, condition, breaks=100, no.zero=TRUE, freq=FALSE, ...)
{
    if (length(COUNT)>0)
    {
        ZERO <- length(COUNT[COUNT == 0])
        if (no.zero) {
            COUNT <- COUNT[COUNT != 0]
        }
        hist(COUNT, breaks, xlab = "intensities", freq=freq, main = paste("Histogram of" , condition), ...)
        rug(COUNT)
        
        COUNT.NA <- COUNT
        COUNT.NA[COUNT.NA == -Inf] <- NA
        mn <- mean(COUNT.NA, na.rm = TRUE)
        stdev <- sd(COUNT.NA, na.rm = TRUE)
        x <- COUNT
        curve(dnorm(x, mean = mn, sd= stdev), add=TRUE, col="red", lty="dotted", xaxt="n")
        quant <- quantile(COUNT.NA, na.rm = TRUE)
        abline(v=quant[2], col="blue")
        abline(v=quant[3], col="red")
        abline(v=quant[4], col="blue")
        
        mtext(paste("mean ", round(mn, 1), "; sd ", round(stdev, 1), "; N ", length(COUNT), "; N zero ", ZERO, sep=""), side=1, cex=.75)
        return(quant)
    } 
} 

## Row wise normalisation
rowshift <- function(df, batch, ref=NULL, useMean = FALSE) {
    
    # df a matrix where rows are genes/proteins and columns are samples
    # batch a factor or vecter has the same length as ncol(df) 
    # to indicate the batch assignment of samples
    
    # ref a logical vector has the same length as ncol(df) to indicated which
    #   columns are the common references among batches. 
    # If it is NULL (by default), the mean of all channels will be used as batch reference.
    feat <- rownames(df)
    
    if (is.data.frame(df)){
        df <- apply(df, 2, as.numeric)
    }
    
    if (is.null(ref)){
        ref <- 1:ncol(df)
    }
    
    df[df == 0] <- NA
    df <- log2(df)
    df[df == -Inf] <- NA
    
    
    
    b_ref <- batch[ref]
    expr_ref <- df[, ref, drop=FALSE]
    grandmeans <- rowMeans(expr_ref, na.rm = TRUE)
    grandmeans2 <- rowMeans(df, na.rm = TRUE)
    
    for (i in unique(batch)) {
        batchmean <- rowMeans(expr_ref[, b_ref == i, drop = FALSE], na.rm = TRUE)
        off <- batchmean - grandmeans
        
        if (useMean) {
            batchmean <- rowMeans(df[, batch == i, drop = FALSE], na.rm = TRUE)
            off2 <- batchmean - grandmeans2
            ona <- is.na(off)
            off[ona] <- off2[ona]
        }
        df[, batch == i] <- df[, batch == i] - off
    }
    df <- as.data.frame(2^df)
    rownames(df) <- feat
    df[is.na(df)] <- 0
    return(df)
}


## Row wise normalisation OLD
rowwise_normalise <- function(data, quantil_cutoff = TRUE, add.jitter = TRUE, shift1 = FALSE) {
    ## Hannes normalization
    # normalize runs (replicates are labeled differently)
    
    # remove rows with 0s only
    data<- data[!rowSums(!data)==ncol(data), ]
    # histnorm3(rowSums(step3.data==0), "data zero blocks distribution", freq = T)
    
    # continue with data for mean calculation with NA instead of 0
    data.NA <- data
    data.NA[data.NA == 0] <- NA 
    
    
    # get position of PS for run normalization
    PS1.pos <- grep("PS1", colnames(data))
    PS2.pos <- grep("PS2", colnames(data))
    PS.pos <- c(PS1.pos, PS2.pos)
    runs <- length(PS1.pos)
    
    # calculate row means
    all_PS_row_means<- rowMeans(data.NA[, PS.pos], na.rm = TRUE, dims = 1)
    # all_PS_row_means[is.na(all_PS_row_means)]
    # histnorm3(log(all_PS_row_means), "all PS row means", freq = T)
    
    # test for na in data
    if (!isempty(data[is.na(data)])) {
        warning("data contains non numeric elements")
    }
    
    data.norm <- data
    
    for(i in 1:runs) {
        j <- i-1
        data.runwise.PS <- cbind(data.NA[ ,PS1.pos[i]], data.NA[ , PS2.pos[i]])
        row.names(data.runwise.PS) <- row.names(data.norm)
        runwise_PS_row_means <- rowMeans(data.runwise.PS, na.rm = TRUE, dims = 1)
        
        # shift data 1 up (remove zeros)
        if (shift1) {
            runwise_PS_row_means[is.na(runwise_PS_row_means)] <- 0
            runwise_PS_row_means <- runwise_PS_row_means + 1
        }
        
        norm_factor <- all_PS_row_means/runwise_PS_row_means
        # histnorm3(norm_factor, paste("normfactor before outlier elimination: run ",  i))
        
        
        # test for na's in dataset
        if (!isempty(norm_factor[is.na(norm_factor)])) {
            warning(paste("norm_factor contains NAs:", length(norm_factor[is.na(norm_factor)])))
        }
        
        
        norm_factor.nacount <- length(norm_factor[is.na(norm_factor)])
        if (!isempty(norm_factor[is.infinite(norm_factor)])) {
            warning(paste("norm_factor contains infinite elements:", length(norm_factor[is.infinite(norm_factor)])))
        }
        
        
        norm_factor.infcount <- length(norm_factor[is.infinite(norm_factor)])
        
        # quantile cutoff
        if (quantil_cutoff) {
            quant <- quantile(norm_factor, c(0, 0.05, 0.50, 0.95), na.rm = TRUE)
            
            norm_factor_mod <- matrix(0, length(norm_factor), 1)
            mod <- which(is.na(norm_factor) | norm_factor < quant[2] | norm_factor > quant[4] | is.infinite(norm_factor))
            norm_factor_mod[mod] <- 1
            
            
            if(add.jitter) {
                quant3.na.jit <- jitter(rep(quant[3], norm_factor.nacount), factor = 3, amount = NULL)
                quant3.inf.jit <- jitter(rep(quant[3], norm_factor.infcount), factor = 3, amount = NULL)
                norm_factor[is.na(norm_factor)] <- quant3.na.jit
                # (infinite values should not appear any more)
                norm_factor[is.infinite(norm_factor)] <- quant3.inf.jit
                quant2.jit <- jitter(rep(quant[2], sum(norm_factor < quant[2])), factor = 3, amount = NULL)
                quant4.jit <- jitter(rep(quant[4], sum(norm_factor > quant[4])), factor = 3, amount = NULL)
                norm_factor[norm_factor < quant[2]] <- quant2.jit
                norm_factor[norm_factor > quant[4]] <- quant4.jit
                
            } else {
                norm_factor[is.na(norm_factor)] <- quant[3]
                norm_factor[norm_factor < quant[2]] <- quant[2]
                norm_factor[norm_factor > quant[4]] <- quant[4]
                norm_factor[is.infinite(norm_factor)] <- quant[3]
            }
        }
        
        else {     # simple normfactor gap filling by 1
            norm_factor[is.na(norm_factor)] <- 1
            norm_factor[is.infinite(norm_factor)] <- 1  
        }
        
        
        if (!isempty(norm_factor[is.na(norm_factor)])) {
            warning("norm_factor contains non numeric elements")
        }
        
        data.norm[, (1+j*10):(10+j*10)] <- apply(data.norm[, (1+j*10):(10+j*10)], 2, "*", norm_factor)
        # histnorm3(norm_factor, paste("normfactor after outlier elimination: run ",  i), breaks = 100)
    }
    
    # in rows with only 0 counts/intensities a division by 0 leads to NaN (should not appear any more)
    if (!isempty(data.norm[is.na(data.norm)])) {
        warning("data.norm contains non numeric elements")
    }
    
    # histnorm3(as.matrix(data.norm), "row wise normalised")
    return(data.norm)
} 





# Between sample normalisation
colwise_normalise <- function(data) {
    # normalize samples (replicates are labeled equal)
    
    # for analysis float values are converted to int (count data)
    maxint<- .Machine$integer.max
    data.scale<-as.data.frame(rescale(as.matrix(data), to = c(0, maxint))) 
    data.int <- round(data.scale, digits=0)
    
    #size factors
    sf= DESeq2::estimateSizeFactorsForMatrix(data.int)
    # barplot(sf, las=2)
    
    # normalized count data
    data.norm <- as.data.frame(t(t(data.int)/sf))
    
    # histnorm3(data.norm, "data DESeq normalised")
    return(data.norm)
} 


## Impute with minimal value across the dataset
impute_min <- function(data, annotation) {
    data.NA <- data
    data.NA[data.NA == 0] <- NA
    
    global_minimum <- min(data.NA, na.rm = T) 
    imputation_value <- round(global_minimum) - 1
    
    data.imp <- data.frame("geneID" = as.character(rownames(data)) )
    data <- as.data.table(data)
    
    annotation <- annotation[SAMPLE_ID %in% colnames(data) ]
    for ( i in unique(annotation$PROTEOMICS_BATCH) ){
        samples <- unique(annotation[PROTEOMICS_BATCH == i]$SAMPLE_ID)
        df <- data[, ..samples ]
        df$protZeroFreq <- apply(df,1,function(x){ sum(x == 0) })
        df[protZeroFreq <= 3 & df == 0] <- imputation_value
        df$protZeroFreq <- NULL
        data.imp <- cbind(data.imp, df)
    }
    data.imp <- as.data.frame(data.imp)
    data.imp$geneID <- as.character(data.imp$geneID)
    return(data.imp)
}


normalize_expression_matrix <-  function(ematrix, sizefactor=FALSE, rowcenter=TRUE, log2scale=TRUE, nonzero=0
){
    require(DESeq2)
    
    if(sizefactor){
        sf= DESeq2::estimateSizeFactorsForMatrix(ematrix)
        ematrix= t(t(ematrix)/sf)
    }
    
    if(log2scale){
        ematrix= log2(nonzero + ematrix)
        ematrix[which(!is.finite(ematrix))] <- NA  
    }
    
    if(rowcenter){
        if(log2scale){
            ematrix= ematrix - rowMeans(ematrix, na.rm=T)
        }else{
            if(any(ematrix <0, na.rm=T)){
                ematrix= ematrix - rowMeans(ematrix, na.rm=T)
            }else{
                ematrix= ematrix / rowMeans(ematrix, na.rm=T)
            }
        }
    }
    return(ematrix)
}



get_proteome_design_matrix_simple <- function(
    patient_id = NULL,
    all_sample_ids = NULL,
    column_name = 'SAMPLE_ID',
    covariates = NULL
){
    design_mat <- matrix(0, 
                         nrow=length(all_sample_ids), 
                         ncol=2, 
                         dimnames=list(all_sample_ids, c('Intercept', column_name))
    )
    design_mat[, "Intercept"]=1
    design_mat[patient_id, column_name]=1
    
    
    if(!is.null(covariates)){
        covariates <- covariates[match(rownames(design_mat), rownames(covariates)),] #rownames should contain ids
        design_mat <- cbind(design_mat, as.matrix(covariates))
    }
    return(design_mat)  
}



get_protein_limma_differential_expression <- function(
    protein_intensity_mat,
    design_matrix,
    coeff_patient = "SAMPLE_ID"
){
    require(limma)
    require(data.table)
    
    #' initial Fit
    lfit= lmFit( protein_intensity_mat, design = design_matrix)
    
    #' * perform moderated t-Test: compute P-values with empirical Bayes
    mod_tstats <- eBayes( contrasts.fit( lfit, coefficients = coeff_patient) )
    
    #' * get fold changes of interest
    prot_limma_de_res= as.data.table(topTable(
        mod_tstats, coef = coeff_patient,
        number=Inf, confint=T, sort.by='none', genelist=rownames(mod_tstats)
    ))
    setnames(prot_limma_de_res, c(
        'GENEID', 'PROTEIN_LOG2FC', 'prot_ci_left', 'prot_ci_right', 'prot_baseMean',
        'prot_mod_t_stat', 'PROTEIN_PVALUE', 'PROT_PADJ', 'prot_logodds' ))
    
    return(prot_limma_de_res)
}


compute_protein_zscore <- function(
    protein_limma_dt,
    column_intensity = 'PROTEIN_LOG2INT',
    column_foldchange = 'PROTEIN_LOG2FC'
){
    
    # compute std deviation per gene
    protein_limma_dt[, SD_PROTEIN_LOG2INT := sd(get(column_intensity), na.rm=T), by= GENEID]
    
    # compute Z-score
    protein_limma_dt[, PROTEIN_ZSCORE := get(column_foldchange)/SD_PROTEIN_LOG2INT]
    return(protein_limma_dt)
}


wrapper_aberrant_protein_expr_simple <- function(
    prot_intensity,
    coln_sample_id = 'SAMPLE_ID',
    p_adjust_method = 'hochberg',  
    LIMMA_covariates = NULL,
    normalize_matrix = TRUE
){
    library(data.table)
    
    # normalize input by size factors
    if(normalize_matrix){ 
        prot_log2_norm_intensity_mat = normalize_expression_matrix(prot_intensity, sizefactor = T, rowcenter = F, log2scale = T, nonzero = 0)
        prot_log2_norm_intensity_mat[which(!is.finite(prot_log2_norm_intensity_mat))] <- NA  
    }else{
        prot_log2_norm_intensity_mat = prot_intensity
        prot_log2_norm_intensity_mat[which(!is.finite(prot_log2_norm_intensity_mat))] <- NA  
    }
    
    
    
    # call limma with specific design for each sample
    res <- sapply(colnames(prot_log2_norm_intensity_mat), 
                  function(sample_id){
                      mydesign <- get_proteome_design_matrix_simple(
                          patient_id = sample_id, 
                          all_sample_ids = colnames(prot_log2_norm_intensity_mat),
                          column_name = paste0(coln_sample_id, 'case'),
                          covariates = LIMMA_covariates
                      )
                      # limma fit
                      single_proteome_limma_res= get_protein_limma_differential_expression(
                          prot_log2_norm_intensity_mat, 
                          design_matrix = mydesign,
                          coeff_patient = paste0(coln_sample_id, 'case')
                      )
                      
                      # multiple testing correction
                      single_proteome_limma_res[,PROTEIN_PADJ := p.adjust(PROTEIN_PVALUE, method = p_adjust_method)]
                      
                      # add sample
                      single_proteome_limma_res[,eval(coln_sample_id):= sample_id]
                      
                      # update res
                      setnames(single_proteome_limma_res, toupper(names(single_proteome_limma_res)))
                      return(single_proteome_limma_res)
                  }, 
                  simplify = FALSE
    )
    resdt <- do.call("rbind", res)

    # get tidy normalized expression
    prot_norm_intensity_dt <- melt(
        as.data.table(prot_log2_norm_intensity_mat, keep.rownames = T), 
        id.vars='rn',
        measure.vars = colnames(prot_log2_norm_intensity_mat),
        value.name = 'PROTEIN_LOG2INT',
        variable.name = toupper(coln_sample_id)
    )
    setnames(prot_norm_intensity_dt,'rn', 'GENEID')
    
    # get rank on normalized intensity
    prot_norm_intensity_dt[which(!is.finite(prot_norm_intensity_dt$PROTEIN_LOG2INT)),]$PROTEIN_LOG2INT <- NA

    
    # merge norm intensity with limma 
    resul <- as.data.table(merge(prot_norm_intensity_dt, resdt))
    # compute Z-score
    resul <- compute_protein_zscore(resul)
    resul <- resul[!duplicated(resul), ]
    resul[, PROTEIN_outlier := PROTEIN_PADJ < .1]
    return(resul)
}


plotAberrantProteinPerSample <- function(ods, main, padjCutoff=0.05, zScoreCutoff=0,
                                         outlierRatio=0.001, col=brewer.pal(3, 'Dark2')[c(1,2)],
                                         yadjust=c(1.2, 1.2), labLine=c(3.5, 3), ymax=NULL,
                                         ylab="#Aberrantly expressed proteins", labCex=par()$cex, ...){
    
    if(missing(main)){
        main <- 'Aberrant Proteins per Sample'
    }
    count_vector <- ods$N
    names(count_vector) <- ods$SAMPLE_ID
    count_vector <- sort(count_vector)
    
    ylim <- c(0.4, max(1, count_vector)*1.1)
    if(!is.null(ymax)){
        ylim[2] <- ymax
    }
    
    replace_zero_unknown <- 0.5
    ticks <- c(replace_zero_unknown, signif(10^seq(
        from=0, to=round(log10(max(1, count_vector))), by=1/3), 1))
    
    labels_for_ticks <- sub(replace_zero_unknown, '0', as.character(ticks))
    
    bp= barplot2(
        replace(count_vector, count_vector==0, replace_zero_unknown),
        log='y', ylim=ylim, names.arg='', xlab='', plot.grid=TRUE,
        grid.col='lightgray', ylab='', yaxt='n', border=NA, xpd=TRUE,
        col=  "gray70",  
        main=main)
    
    n_names <- floor(length(count_vector)/20)
    xnames= seq_len(length(count_vector)) #seq_len(n_names*20)
    axis(side=1, at= c(0,bp[xnames,]), labels= c(0,xnames))
    axis(side=2, at=ticks, labels= labels_for_ticks, ylog=TRUE, las=2)
    
    # labels
    mtext('Sample rank', side=1, line=labLine[1], cex=labCex)
    mtext(ylab, side=2, line=labLine[2], cex=labCex)
    
    # legend and lines
    hlines = c(Median=ifelse(median(count_vector)==0, replace_zero_unknown,
                             median(count_vector)) , Quantile90=quantile(
                                 count_vector, 0.9, names=FALSE))
    color_hline= c('red','black')
    abline(h= hlines, col=color_hline)
    text(x=c(1,1), y= hlines*yadjust, col=color_hline, adj=0,
         labels=c('Median', expression(90^th ~ 'percentile')))
    
    box()
}


replace_with_mean <- function(vec) {
    m <- mean(vec, na.rm = TRUE)
    vec[is.na(vec)] <- m
    return(vec)
}






# Plots for normalisation progress

plots_normalization <- function(data, annotation, remove_samples = NULL) {
    #________________
    ## prepare data for heatmap-creation - raw intensities all proteins without zeroblocks
    
    # barplots
    dfx <- data.table( SAMPLE_ID = rownames(t(data)) , SUM_INT = rowSums(t(data)) )
    dfx <- merge(annotation[ ,.(SAMPLE_ID, PROTEOMICS_BATCH, BATCH_RUN, INSTRUMENT)], dfx, by = "SAMPLE_ID" )
    dfx <- dfx[order(BATCH_RUN)]
    
    ggplot(dfx, aes(SAMPLE_ID , SUM_INT))+
        geom_col()+ 
        theme_bw()+ 
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
              axis.title.x = element_blank())+
        facet_wrap(~PROTEOMICS_BATCH , scales = "free_x", nrow =2 )
    
    
    # Heatmap and PCA
    data <- data[, !(names(data) %in% remove_samples)]
    data[is.na(data)] <- 0
    data <- log2(as.matrix(data))
    data[which(!is.finite(data))] <- 0 
    
    annotation <- as.data.frame(annotation)
    rownames(annotation) <- annotation$SAMPLE_ID
    annotation$INSTRUMENT <- as.character(annotation$INSTRUMENT)
    pheatmap(cor(data), annotation_col = annotation[ , c("PROTEOMICS_BATCH", "BATCH_RUN",  "INSTRUMENT")])
    
    
    
    ## PCA of not normalized data
    data <- data[ which(apply(data, 1, var) !=0 ), ]
    pca <-prcomp(t(data), scale = T, center = T)
    pca.var<-pca$sdev^2/sum(pca$sdev^2)
    pca_df <- data.frame(x = pca$x[,1], y = pca$x[,2])
    pca_df <- merge(pca_df, annotation,  by=0)
    
    ggplot(pca_df, aes(x, y, color = PROTEOMICS_BATCH, shape =  INSTRUMENT )) +
        geom_point() +
        theme_bw() +
        theme(aspect.ratio=1) 
}







