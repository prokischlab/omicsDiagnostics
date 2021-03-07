############################################
### autoencoder models
############################################

pacman::p_load(data.table, keras, SummarizedExperiment, pcaMethods )


### special wraper for SummarizedExperiment object with shuffle
run_autoencoder_model_wrapper = function(model_type, se, ...) {
    shuf_se = se[,sample(ncol(se))]
    
    ae = run_autoencoder_model(model_type = model_type, X = t(assay(shuf_se, 'X')), X_corr = t(assay(shuf_se, "X_corr")),
                                X_na = t(assay(shuf_se, "X_na")), batches = colData(shuf_se), ...)
    
    sortback = match(colnames(se), colnames(shuf_se))
    pred = t(ae$predicted[sortback,])
    
    assays(se)$X_pred = pred
    assays(se)$X_pred[assays(se)$X_na] = 0
    
    return(list(predicted = pred, model = ae$model, other=ae$other, se=se))
}  


### wrapper for different ae models 
run_autoencoder_model = function(model_type, X, X_corr, X_na, batches, encoding_dim, seed, epochs, batch_size, es_patience, optimizer, verbose, learning_rate) {
  ae_model = switch(model_type,
                     "AE_batches" = run_autoencoder(X, X_corr, NA_mask = X_na, encoding_dim = encoding_dim, batch_cov = batches, seed = seed, epochs=epochs, batch_size = batch_size, patience=es_patience, optimizer=optimizer, verbose = verbose, learning_rate = learning_rate),
                     "AE_batches_Dfit" = run_autoencoder_Dfit(X, X_corr, NA_mask = X_na, encoding_dim = encoding_dim, batch_cov = batches, seed = seed, epochs=epochs, batch_size = batch_size, patience=es_patience, optimizer=optimizer, verbose = verbose, learning_rate = learning_rate), 
                     "AE_woBatches" = run_autoencoder(X, X_corr, NA_mask = X_na, encoding_dim = encoding_dim, batch_cov = NULL, seed = seed, epochs=epochs, batch_size = batch_size, patience=es_patience, optimizer=optimizer, verbose = verbose, learning_rate = learning_rate),
                     
                     ### pca with batches does not make sense
                     "PCA" = run_PCA(X, pcs=encoding_dim, batch_cov=NULL, NA_mask=X_na),
  )
  return(list(predicted = ae_model$predicted, model = ae_model$model, other=ae_model))
}






# 
# X = t(assays(se)$X)
# X_corr = t(assays(se)$X_corr)
# NA_mask = t(assays(se)$X_na)
# batch_cov = colData(se)
# encoding_dim = 22
# seed=NULL
# epochs=100
# batch_size=500
# optimizer="Adam"
# learning_rate = 0.001
# patience=500
# verbose=1









# run pca analysis on protein intensities and reports predicted values
run_PCA = function(inten, pcs=10, batch_cov=NULL, NA_mask=NULL) {
  if(is.null(NA_mask)){
    inten = na.omit(inten)
  } else {
    inten[NA_mask] = NA
  }
  pca_matrix = inten
  
  if(!is.null(batch_cov)) {
    pca_matrix = cbind(pca_matrix, batch_cov)
  }
  
  # actual pca
  pca_prot = pca(pca_matrix, nPcs=pcs)
  pca_pred = predict(pca_prot, pca_matrix)$x
  
  if(!is.null(batch_cov)) {
    pca_pred = pca_pred[,1:(ncol(pca_pred)-ncol(batch_cov))] # removed batch columns at end
  }
  
  # remove missing values
  if(!is.null(NA_mask)) { 
    pca_pred[is.na(pca_matrix[,1:ncol(inten)])] = 0 }
  
  dimnames(pca_pred) = dimnames(inten)
  return( list(predicted = pca_pred, model = pca_prot))
}





# Run covariate-aware keras autoencoder with one hidden layer to get expected intensities.
run_autoencoder = function(inten, inten_corr, NA_mask, encoding_dim, batch_cov=NULL,
                           seed=NULL, epochs=300, batch_size=8, optimizer="Adam", learning_rate = 0.001,
                           patience=10, verbose=1) {
  
  NA_mask = as.matrix(NA_mask)
  inten = as.matrix(inten)
  inten_corr = as.matrix(inten_corr)
  
  ### normally not necessary but want to make sure
  inten[NA_mask] = 0
  inten_corr[NA_mask] = 0
  
  if(!is.null(seed)){
    set.seed(seed)
    tensorflow::use_session_with_seed(seed)
  }
  
  ## proteins
  p = ncol(inten)
  
  ### removed split to not further decrease sample size
  # ### Train/test split
  # train = sample(nrow(inten), floor(nrow(inten)*0.8))
  # x_corr_train = inten_corr[train,]
  # y_train = inten[train,]
  # x_corr_test = inten_corr[-train,]
  # y_test = inten[-train,]
  
  # train = sample(nrow(inten))
  x_corr_train = inten_corr#[train,]
  y_train = inten#[train,]
  
  if(!is.null(batch_cov)) {
    batch_cov = as.matrix(batch_cov)
    batch_train = batch_cov#[train,]    
    #batch_test = batch_cov[-train,]
    input_fit = list(cbind(x_corr_train, batch_train), batch_train)
    # input_evaluate = list(cbind(x_corr_test,batch_test), batch_test)
    input_predict = list(cbind(inten, batch_cov), batch_cov)
  } else {
    input_fit = x_corr_train
    # input_evaluate = x_corr_test
    input_predict = inten
  }
  
  ######### MODEL ##########
  # create and compile model
  base_model = get_keras_model(inten, batch_cov, encoding_dim)
  
  es = callback_early_stopping(monitor = "loss", patience = patience,
                                restore_best_weights = TRUE, min_delta = 0.001)
  
  opti = switch(optimizer,
                 "Adam" = optimizer_adam(lr = learning_rate),
                 "rmsprop" = optimizer_rmsprop(lr = learning_rate),
                 "sgd"=optimizer_sgd(lr=learning_rate, nesterov=TRUE))
  
  base_model %>% compile(
    #optimizer = 'rmsprop',
    optimizer = opti,
    loss = 'mean_absolute_error',
    # loss = 'mean_squared_error',
    metrics = c('mean_squared_error')
  )
  # ### Fit the model on whole data
  history = base_model %>% keras::fit(
    x = input_fit,
    y = y_train,
    verbose=verbose,
    epochs=epochs, batch_size=batch_size,  
    callback = es
  )
  
  # plot(history, main = "History of fitting the model, batch-aware PCA", method="ggplot2") + scale_y_log10()
  #eval = base_model %>% evaluate(input_evaluate, y_test) ### evaluation not necessary
  
  
  X_hat = base_model %>% predict(input_predict)
  dimnames(X_hat) = dimnames(inten)
  X_hat[NA_mask] = 0
  
  list(predicted = X_hat, model=base_model, epoch=es$stopped_epoch, history=history)
}



### returns appropriate autoencoder model with/without batch effects considered
### input (+BE) -> latent_pca (+BE) -> output
get_keras_model = function(x, batches, encoding_dim, encod_use_bias=FALSE) {
  if(is.null(batches)) {
    input = layer_input(shape = ncol(x))
    masking_NAs = input %>% layer_masking(mask_value = 0)
    # latent space of pca
    latent_pca = masking_NAs %>% layer_dense(units = encoding_dim, use_bias=encod_use_bias)
    output = latent_pca %>% layer_dense(units = ncol(x))
    base_model = keras_model(inputs = input, outputs = output) 
  } else {
    # input layer = protein + batch 
    input = layer_input(shape = c(ncol(x)+ncol(batches)))
    # NAs have been substituted for 0 
    masking_NAs = input %>% layer_masking(mask_value = 0)
    # latent space of pca
    latent_pca = masking_NAs %>% layer_dense(units = encoding_dim, use_bias=encod_use_bias, activation = "linear" ) 
    cov_input = layer_input(shape = ncol(batches))
    output = layer_concatenate(c(latent_pca, cov_input))  %>%  layer_dense(units = ncol(x), activation = "linear" )
    base_model = keras_model( inputs = c(input, cov_input), outputs = output ) 
  }
  return(base_model)
}




### after fitting: last refinement of decoder per protein 
run_autoencoder_Dfit = function(X, X_corr, NA_mask, encoding_dim, batch_cov, verbose, ...) {
  ae_obj = run_autoencoder(X, X_corr, NA_mask = NA_mask, encoding_dim = encoding_dim, batch_cov = batch_cov, verbose = verbose, ... )
  
if(verbose) {
    print('before decoder refinement: ')
    x_pred_ae = ae_obj$predicted
    x_pred_ae[NA_mask] = NA
    print(paste0('mae: ', mean_absolute_error(x_pred_ae, X) ))
  }
  
  ae_model = ae_obj$model
  weights = get_weights(ae_model)
  
  B = as.matrix(batch_cov)  # batch_cov : samples x con
  H = as.matrix(cbind(X,B)) %*% weights[[1]]   # weights[[1]] : prot x q
  cbi = as.matrix(cbind(H,B)) # cbi : samples x (q+con)
  
  ### check linear model refinement fit 
  # x_decod = X
  # x_pred = as.data.frame(ae_obj$predicted)
  # x_pred[NA_mask] = NA
  # # which.max(apply(decod_weights, 2, max, na.rm=TRUE))
  # par(mfrow=c(1,2))
  # for( i in c(4, 50, 5141, 3405, 6843)) {#,1378,5141)) {
  #     prot_num = i # 10
  #     prot_coef = coef(lm(x_decod[, prot_num]~cbi) )
  #     prot1 = lm(x_decod[, prot_num] ~ cbi)
  #     plot(na.omit(x_decod[,prot_num]), predict(prot1), main=colnames(x_decod)[prot_num], ylim=c(-3,6), xlim=c(-3,6))
  #     abline(0,1,col="blue")
  #     plot(x_decod[,prot_num], x_pred[,prot_num], main="without lm", ylim=c(-3,6), xlim=c(-3,6))
  #     abline(0,1,col="blue")
  # }
  # par(mfrow=c(1,1))
  
  
  ### linear regression for each protein, mask NA
  X[NA_mask] = NA
  decod_matrix = apply(X, 2, function(x) {
    x_pos = which(!is.na(x))
    coef(lm(x[x_pos]~cbi[x_pos,]) )
  })
  
  ### set new decoding weights
  decod_bias = decod_matrix[1,]
  dim(decod_bias) = length(decod_bias)  # important to get into weights
  decod_weights = decod_matrix[-1,]
  decod_weights[is.na(decod_weights)] = 0  # na weights to 0
  
  new_weights = list(weights[[1]], unname(decod_weights), unname(decod_bias) )
  set_weights(ae_model, weights = new_weights)  # set new weights to ae model
  
  ### updated prediction with new weights
  X_pred_input = X
  X_pred_input[NA_mask] = 0 
  input_predict = list(cbind(X_pred_input, B), B)
  X_pred = ae_model %>% predict(input_predict)
  dimnames(X_pred) = dimnames(X)
  X_pred[NA_mask] = NA
  
  if(verbose) {
    print('after decoder refinement: ')
    print(paste0('mae: ', mean_absolute_error(X_pred, X) ))
  }
  
  return ( list(predicted = X_pred, model=ae_model) )
  
}






