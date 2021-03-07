############################################
### data input methods
############################################


summarize_prot = function(file_proteins, file_annotation, confounders_used){
  ### protein intensity data, samples x genes
  X_raw = as.matrix(t(file_proteins))
  
  ### get batches to one-hot-encoding
  batch_table = input_batches_to_oneh(file_annotation, confounders_used = confounders_used)
  
  ### sort samples
  input_batches_order = match(rownames(X_raw), rownames(batch_table))
  batch_oneh_table = batch_table[input_batches_order,]
  
  ### get NA data and normalise
  X_na = X_raw == 0
  
  sf <- DESeq2::estimateSizeFactorsForMatrix(t(X_raw))
  X_raw_trans <- X_raw/sf
  X <- log2(X_raw_trans)
  X[X_na] = 0
  X_c = center_columnwise(X, X_na)
  
  se = SummarizedExperiment(assays = list(X=t(X_c), X_raw = t(X_raw), X_na = t(X_na), X_mean=t(X-X_c)),
                            colData = batch_oneh_table)
  metadata(se)$batches = file_annotation
  metadata(se)$size_factors = sf
  return(se)
}





### center data to zero per protein, for norm data
center_columnwise = function(X_raw, X_na) {
  X_norm = X_raw
  X_norm[X_na] = NA
  X_norm = t(t(X_norm) - colMeans(X_norm, na.rm = T))
  X_norm[X_na] = 0
  return(X_norm)
}



input_batches_to_oneh = function(input_batches, confounders_used){
  inp = data.frame(input_batches)
  #rownames(inp) = inp$SAMPLE_ID
  inp = inp[,confounders_used]
  
  one_hot_list = lapply(confounders_used, function(x) {
    prok_batch = factor(inp[[x]])
    prok_batch_cov = to_categorical(as.integer(prok_batch)-1, num_classes = nlevels(prok_batch))  ## error source
    ### error if tensorflow/keras is not correctly installed, but then training also does not work
    colnames(prok_batch_cov) = paste0(x,"_",levels(prok_batch))
    prok_batch_cov
  })
  
  outDf = do.call(cbind, one_hot_list)
  rownames(outDf) = rownames(inp)
  colnames(outDf) = tolower( colnames(outDf))
  return(outDf)
}

















############################################
#' create some simulated data
get_simulated_data = function( num_prot=10000, num_batches = 14, num_samples_per_batch = 10) {
  
  p = num_prot  # proteins
  b  = num_batches  # batches
  n.per.b = num_samples_per_batch  # samples per bacth
  n = b*n.per.b
  
  batch = factor(rep(paste0("batch_", 1:b), each = n.per.b))    # batches 
  batch_oneh_sim = to_categorical(as.integer(batch)-1, num_classes = nlevels(batch))   # one-hot encoding of batches
  
  ### simulate batch effects 
  
  beta_b = matrix(rnorm(p*b), nc=b)   # (protein x batches)
  B = t(beta_b[,as.integer(batch)])   # (samples x proteins)
  
  # q = true size of the latent space
  q = 5
  H = matrix(rnorm(n*q), nc=q)
  D = matrix(rnorm(q*p), nc=p)   # decoder
  
  ## Simulate noise
  # standard dev (=1 for each protein)
  sd = rep(1,p)
  eps = matrix(sd*rnorm(n*p),nc=p,byrow = FALSE)
  
  ## Simulate expression values
  X_sim = H%*%D + B + eps
  
  ### transform to matrix
  X_sim = as.matrix(X_sim)
  dimnames(X_sim) = list(paste0("sample_",1:nrow(X_sim)), paste0("prot_",1:ncol(X_sim)))
  
  batch_oneh_sim = as.matrix(batch_oneh_sim)
  dimnames(batch_oneh_sim) = list(paste0("sample_",1:nrow(batch_oneh_sim)), paste0("batch_",1:ncol(batch_oneh_sim)))
  
  X_na = X_sim == 0  # should not occur 
  
  ### samples in columns and colData, proteins in rows !!!
  # https://bioconductor.org/packages/release/bioc/vignettes/SummarizedExperiment/inst/doc/SummarizedExperiment.html
  se = SummarizedExperiment(assays=list(X = t(X_sim), X_na = t(X_na)), colData=batch_oneh_sim)
  return(se)
  #return(list(X = X_sim, batches = batch_oneh_sim))
}



############################################
#' read in old data 
get_prok_raw = function(file_proteins, file_annotation, confounders_used, with_control = FALSE, colname_filter_sample = "USE_FOR_PROTEOMICS_PAPER"){
  input_annotation = data.frame(fread(file_annotation))
  samples_to_keep = input_annotation[input_annotation[, colname_filter_sample] == TRUE, ]$proteome_ID
  
  ### protein intensity data
  raw_data = data.frame(fread(file_proteins), row.names = 1) # TODO EDIT
  
  # with | without controls
  if(with_control) {
    samples_to_keep = c(samples_to_keep, colnames(raw_data)[grepl("PS", colnames(raw_data))] )
  }
  
  X_raw = as.matrix(t(raw_data[,colnames(raw_data) %in% samples_to_keep])) # only keep appropriate samples, samples x genes !
  
  ### remove proteins with too many NaNs
  MIN_SAMPLE_NUM_PROT = 0.7
  prots_to_keep  = colSums(X_raw!=0) >= nrow(X_raw) * MIN_SAMPLE_NUM_PROT 
  
  print(paste0('proteins in >',nrow(X_raw) * MIN_SAMPLE_NUM_PROT,' samples (', MIN_SAMPLE_NUM_PROT,'%):') )
  print(table(prots_to_keep))
  
  X_raw = X_raw[,prots_to_keep]
  
  ### get batches to one-hot-encoding
  input_batches = input_annotation[input_annotation[,colname_filter_sample] == TRUE,]
  batch_table = input_batches_to_oneh(input_batches, confounders_used = confounders_used)
  
  ### sort samples
  input_batches_order = match(rownames(X_raw), rownames(batch_table))
  batch_oneh_table = batch_table[input_batches_order,]
  
  
  ### remove full batch if #samples <3
  print( paste0('before batch-specific removal: intensities == NA: ', sum(X_raw==0), ' -> ', round(sum(X_raw==0)/ (dim(X_raw)[1]*dim(X_raw)[2]) ,3) ))
  batch_subset = batch_oneh_table[, grep("proteomics_batch_", colnames(batch_oneh_table))]
  # iterate trough each protein and check number of samples per batch
  X_rmvd = sapply(1:ncol(X_raw), function(i) {
    x = X_raw[,i]
    sa = batch_subset[names(x)[!x==0],, drop=FALSE]   # get samples per prot
    samples_per_mix = apply(sa,2, sum)
    mix_del = samples_per_mix < 3 #min 3 per batch
    
    m_del = batch_subset[,which(mix_del),drop=FALSE]  # if #mix ==1
    m_del_num = names(which(rowSums(m_del) == 1))
    x[m_del_num] = 0
    x
  } )
  print( paste0('after batch-specific removal: intensities == NA: ', sum(X_rmvd==0), ' -> ', round(sum(X_rmvd==0)/ (dim(X_rmvd)[1]*dim(X_rmvd)[2]) ,3) ))
  dimnames(X_rmvd) = dimnames(X_raw)
  X_raw = X_rmvd
  
  
  ### get NA data and normalise
  X_na = X_raw == 0
  
  sf <- DESeq2::estimateSizeFactorsForMatrix(t(X_raw))
  X_raw_trans <- X_raw/sf
  X <- log2(X_raw_trans)
  X[X_na] = 0
  X_c = center_columnwise(X, X_na)
  
  se = SummarizedExperiment(assays = list(X=t(X_c), X_raw = t(X_raw), X_na = t(X_na), X_mean=t(X-X_c)),
                            colData = batch_oneh_table)
  metadata(se)$batches = input_annotation
  metadata(se)$size_factors = sf
  return(se)
}





