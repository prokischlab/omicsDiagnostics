############################################
### methods to insert noise and artificial outlier to the data
############################################

### corrupt intensities with gaussian noise 
inject_noise = function(se, noise_factor, show_plot=TRUE) {
  inten = t(assays(se)$X)
  NA_mask = t(assays(se)$X_na)
  
  inten = as.matrix(inten)
  inten[NA_mask] = NA

  ### check again - distribution addition seems kind of skrewed TODO
  inten_corr = inten + t( matrix(rnorm(nrow(inten)*ncol(inten)), nc = nrow(inten)) * noise_factor * colSds(inten, na.rm = T))
  inten_corr[is.na(inten)] = 0
  if (show_plot) {
    plot(x =inten, y = inten_corr, main=paste("noise injection to intensities"),
         xlab = "raw", ylab = "corrupted",pch=16)
    grid()
    abline(0,1, color="orange")
  }
  
  assays(se)$X_corr =  t(inten_corr)
  return(se)
}




### inject strong outliers
inject_outliers = function(se, inj_freq = 1e-3, inj_meanlog = 3, inj_sdlog = 1.6, show_plot=TRUE, inj_scheme="both") {

  if("X_woOutlier" %in% names(assays(se)) ) {  # check if outliers already injected
    assays(se)$X = assays(se)$X_woOutlier  
  } 
  
  inten = t(assays(se)$X)
  NA_mask = t(assays(se)$X_na)
  
  inj_direction = switch( inj_scheme,
    "both" = c(0,-1, 1),
    "lower" = c(0,-1, -1),
    "higher" = c(0, 1, 1)
  )
  
  outliersMask = matrix(ncol=ncol(inten), nrow=nrow(inten),
                       sample(inj_direction, prod(dim(inten)), replace=TRUE, 
                              prob=c(1-inj_freq, inj_freq/2, inj_freq/2)))
  dimnames(outliersMask) = dimnames(inten)
  
  inj_zscores = matrix(ncol=ncol(inten), nrow=nrow(inten), 
                       rlnorm(prod(dim(inten)), meanlog=log(inj_meanlog), sdlog=log(inj_sdlog)))
  
  inten[NA_mask] = NA
  inj_zscores[NA_mask] = NA
  
  
  #
  # z * sigma + X_mean = X
  # ==> up/down * z * sigma + X => corruptedIntensity
  # 
  fc = t(t(inj_zscores) * colSds(inten, na.rm=TRUE))
  #print(dim(inten) == dim(outliersMask))
  injections = outliersMask * fc + inten
  
  #
  # ## Check outlier injection
  if(show_plot) {
    table(outliersMask)
    hist(inj_zscores[outliersMask!=0], breaks=100)      
    LSD::heatscatter(inten[outliersMask == -1], injections[outliersMask == -1], 
                     main="Down regulation injections")
    grid()
    abline(0,1)
    
    LSD::heatscatter(inten[outliersMask == 1], injections[outliersMask == 1],
                     main="Upregulation injections")
    grid()
    abline(0,1)
  }
  
  inten[outliersMask != 0] = injections[outliersMask != 0]
  inten[is.na(inten)] = 0
  
  if(show_plot){
    for(i in which(colAnys(outliersMask != 0))[1:4]){
      plot(inten[,i], main=colnames(inten)[i] )
      points(which(outliersMask[,i] != 0), inten[outliersMask[,i] != 0,i], col="red", pch=16)
    }
  }
  
  assays(se)$X_woOutlier = assays(se)$X
  assays(se)$X = t(inten)
  assays(se)$X_out_pos = t(outliersMask)
  
  return(se)
}





