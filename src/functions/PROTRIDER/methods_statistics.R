############################################
### statistics with p-val calculation etc
############################################


### calculate pvalues with different standard deviation approaches
### if with cross validation: runs becomes k-fold
calculate_pvalues = function(se, padj_method = "BY") {  
  observed = as.matrix(t(assays(se)$X) )
  predicted = as.matrix(t(assays(se)$X_pred) )
  NA_mask = t(assays(se)$X_na)
  
  residuals = observed - predicted
  
  if(!is.null(NA_mask)) {
    observed[NA_mask] = NA
    predicted[NA_mask] = NA
    residuals[NA_mask] = NA
  }
  
  # get normal standard deviation with outliers or over multiple runs/cross-validation
  pvalues_sd = colSds( residuals, na.rm = T)
  pvalues = t(pnorm(t(observed), mean = t(predicted), sd = pvalues_sd ))   # assume gaussian
  
  # convert to two-tailed p-values
  pvalues = 2*pmin(pvalues, 1-pvalues)
  dimnames(pvalues) = dimnames(observed)
  
  # multiple testing correction
  applyMargin = 1 # per sample
  pvalues_adj = t( apply( pvalues, applyMargin, p.adjust, method = padj_method))
  
  # z-score calculation
  zscores = t( (t(residuals) - colMeans(residuals, na.rm = T))  / colSds(residuals, na.rm = T) )
  
  assays(se)$X_pval = t(pvalues)
  assays(se)$X_pval_adj = t(pvalues_adj)
  assays(se)$X_z = t(zscores)
  return(se)
}



call_sig_outlier = function(se, fold_change=1, pval_adj=0.1) {
    assays(se)$X_log2fc = get_fold_change(assays(se)$X, assays(se)$X_pred)
    
    assays(se)$X_out_called = (assays(se)$X_pval_adj < pval_adj) &
        (abs(assays(se)$X_log2fc) > log2(fold_change))
    return(se)
}



get_fold_change_raw = function(X, X_pred) {
    fc = X_pred / X
    return(fc)
}

### in log space
get_fold_change = function(X, X_pred) {
    fc = X - X_pred
    return(fc)
}




get_aberrant = function(se, padjCutoff=0.05, by=c("none", "sample", "gene")){
    aberrantEvents = assays(se)$X_pval_adj < padjCutoff
    return(switch(match.arg(by),
                  none = aberrantEvents,
                  sample = colSums(aberrantEvents, na.rm=TRUE),
                  gene = rowSums(aberrantEvents, na.rm=TRUE)
    ))
}






### creates ROC curve and calculates AUC
get_AUC_ROC = function(pValue, isOutlier, show_plot=FALSE, prec_recall = TRUE) {
  score_data = data.table(score=pValue %>% as.vector, 
                           label=(isOutlier!=0) %>% as.vector, 
                           name = "PROTRIDER")
  score_data = na.omit(score_data)  ## TODO: maybe add bootstrapping approach like bootstrapPR()
  pred = prediction( score_data$score, !score_data$label)
  
  if(prec_recall){
    perf = performance(pred,"prec","rec")  # precision prec, recall rec
  } else {
    perf = performance(pred,"tpr","fpr")  # precision prec, recall rec
  }
  
  auc_ROCR = performance(pred, measure = "auc")
  auc_ROCR = auc_ROCR@y.values[[1]]  # AUC value
  
  if(show_plot) {
    plot(perf, main = paste0("AUC: ",round(auc_ROCR,4)), ylim=c(0,1), xlim=c(0,1))
    
      if(prec_recall) {
        abline(1,-1, col="blue")
      } else {
        abline(0,1, col="blue")
      }
  }

    return(list(auc = auc_ROCR, perf_obj = perf))
}






get_prec_rec_auc = function(se){
    scores = -as.vector(assay(se, 'X_pval'))
    labels = as.vector(assay(se, 'X_out_pos') != 0) + 0
    return(get_prec_rec_auc_vector(labels, scores))
}

get_prec_rec_auc_vector = function(labels, scores) {
    library(PRROC)
    labels = labels[!is.na(scores)]
    scores = scores[!is.na(scores)]
    
    pr = pr.curve(scores, weights.class0=labels, curve=TRUE)
    return( list(pr_obj=pr, pr_auc = max(0, pr$auc.integral, na.rm=TRUE)) )
}



root_mean_squared_error = function(fitted, observed){
    sqrt(mean((fitted - observed)^2, na.rm=TRUE))
}

mean_absolute_error = function(fitted, observed){
    mean(abs(fitted - observed), na.rm=TRUE)
}








