############################################
### data handling methods
############################################

### small preview of matrix
pre = function(df) {
    rows = min(nrow(df),6)
    cols = min(ncol(df),6)
    df[1:rows,1:cols]
}


### get a long melted version of one hot encoded batches table - only for one confounding effect
get_batches_melted = function(batches) {
    # batches = batches[, grep("proteomics_batch_", colnames(batches))]
    batches = batches[, grep("batch_", colnames(batches), ignore.case = TRUE)]
    batch_oneh_heatmap = as.data.frame(batches)
    batch_oneh_heatmap$sample = rownames(batches)
    batch_oneh_heatmap$order = 1:nrow(batches)  # problems with rownames ordering 
    batch_oneh_heatmap$reference_sample = NULL
    
    batch_long = reshape2::melt(batch_oneh_heatmap, id = c("sample","order"))
    batch_long = subset(batch_long, value==1)
    batch_long = batch_long[order(as.numeric(batch_long$order)),] 
    rownames(batch_long) = batch_long$sample
    batch_long[ ,c('sample', 'order','value')] = list(NULL)
    colnames(batch_long) = "batches"
    return(batch_long)
}

### melts summarizedExperiment object into long table
### sample | protein | pval | adj pval | sample_protein
get_prot_sample_list = function(se, pval_name ="X_pval", pval_adj_name="X_pval_adj"){
    
    df_pval = assay(se, pval_name)
    long_pval = reshape2::melt(df_pval, id=colnames(df_pval), measure.vars=rownames(df_pval))  # value.name does not work

    outDf = do.call("cbind", list(long_pval, as.vector(assay(se, pval_adj_name)), 
                          as.vector(assay(se, "X_log2fc")),
                          2**as.vector(assay(se, "X_log2fc")), 
                          as.vector(assay(se, "X_out_called")),
                          as.vector(assay(se, "X_raw")),
                          as.vector(assay(se, "X")) + as.vector(assay(se, "X_mean")),
                          as.vector(assay(se, "X")) + as.vector(assay(se, "X_mean")) - as.vector(assay(se, "X_pred")),
                          as.vector(assay(se, "X_z"))
                          ))
    colnames(outDf) = c('protein_id','sample_id','pvalue','adj_pvalue','log2fc','fc','is_outlier', 
                         'intensity_raw','intensity_raw_log2','intensity_norm_log2','z_score')
    

    if(!is.null( assays(se)$X_out_pos )){
        outDf = cbind(outDf, as.vector(assay(se,"X_out_pos")))
        colnames(ncol(outDf)) = "outliers"
    }
        
    outDf = outDf[, !duplicated(colnames(outDf))]
    outDf = na.omit(outDf)
    outDf$sample_id = as.character(outDf$sample_id)
    outDf$protein_id = as.character(outDf$protein_id)
    outDf$sample_prot = paste0(outDf$sample_id, '_',outDf$protein_id)
    return(outDf)
}




transpose_with_names = function(table) {
    oldRN = rownames(table)
    oldCN = colnames(table)
    
    table = t(table)
    rownames(table) = oldCN
    colnames(table) = oldRN
    return(table)
}

















