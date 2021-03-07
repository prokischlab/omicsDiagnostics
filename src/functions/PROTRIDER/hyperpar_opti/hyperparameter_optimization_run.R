############################################
### hyperparamater optimization complete
###########################################


library(reticulate)
use_condaenv('loipfi_protrider')   # conda environment with tensorflwo 1.14.0

source("/data/nasif12/home_if12/loipfins/gitlab/protrider/paper/config.R")  # solve other way

source(file.path(SCRIPT_DIR, "data_input.R") )
source(file.path(SCRIPT_DIR, "data_handling.R") )
source(file.path(SCRIPT_DIR, "methods_noise_outlier.R") )
source(file.path(SCRIPT_DIR, "methods_statistics.R") )
source(file.path(SCRIPT_DIR, "autoencoder_models.R") )

library(ggplot2)


############################################
### get raw data
se = get_prok_raw(DATA_DIR,
                   with_control = WITH_CONTROL_SAMPLES,
                   confounders_used = CONFOUNDERS_USED)
set.seed(SEED)


############################################
### insert noise and artificial outliers
se = inject_outliers(se, inj_meanlog = 3, show_plot=F, inj_scheme="both")
se = inject_noise(se, NOISE_FACTOR, show_plot=F )



############################################
### write autoencoder run with performance into file
check_performance = function(output_file, modelType, se, encoding_dim, batch_size=8, epochs = 300, es_patience=10, seed=NULL, optimizer="Adam", learning_rate=0.001, verbose=0) {
  
  tryCatch({
      ### get model with predicted intensities
      ae_obj = run_autoencoder_model_wrapper(modelType, se, encoding_dim = encoding_dim, seed = seed,
                                              epochs = epochs, batch_size=batch_size, es_patience=es_patience, 
                                              optimizer=optimizer, verbose=verbose, learning_rate=learning_rate)
      se = ae_obj$se
      se = calculate_pvalues(se)
      prec_rec = get_prec_rec_auc(se)
      print(paste0("pre rec auc: ", round(prec_rec[['pr_auc']],5)))

    ### ae history
    if(startsWith(modelType, "PCA")) {
        epochs_stopped = 0
        loss = 0
        val_loss = 0
    } else {
        h = ae_obj$model$history$history
        epochs_stopped = length(h$loss)
        loss = h$loss[[epochs_stopped]]
        val_loss = h$val_loss[[epochs_stopped]]
    }
    
    write(paste0(modelType,",", prec_rec[['pr_auc']],",", loss,",",val_loss,",",epochs_stopped,",",encoding_dim,",",batch_size,",",optimizer,",",learning_rate,",",date()), output_file, append=TRUE)
    
      }, error = function(e) {
        write(paste0(modelType,",", "ERROR",",", "",",","",",","",",",encoding_dim,",",batch_size,",",optimizer,",",learning_rate,",",date()), output_file, append=TRUE)
    } )
}




run_over_param = function(output_file, parallel_cores=4){
    
    ### run only over encoding_dim as hyperparameter
    param_dim = seq(4,66,2)
    param_learningRate = LEARNING_RATE
    param_modelTypes = c("AE_batches_Dfit")
    param_batchSize = BATCH_SIZE
    param_epochs = EPOCHS
    
    ### if PCA was applied
    #all_param = expand.grid(modelType = param_modelTypes, noiseFactor=param_noise, encoding_dim=param_dim, batch_size= param_batchSize, learning_rate = param_learningRate)
    # pca_all_param = data.frame(modelType=c("PCA_batches", "PCA_woBatches"), noiseFactor=0, encoding_dim=param_dim,batch_size=0, learning_rate=0)
    # pca_all_param = data.frame(modelType=c("PCA_batches", "PCA_woBatches"), encoding_dim=param_dim, batch_size=0, learning_rate=0)
    # all_param = rbind(all_param, pca_all_param)  # 10218 x 5 ...
    
    all_param = expand.grid(modelType = param_modelTypes, encoding_dim=param_dim, batch_size= param_batchSize, learning_rate = param_learningRate, epochs=param_epochs)
    

    ### start writing in document
    write("model, auc, loss, val_loss, epochs_stopped, encoding_dim_PC, batch_size, optimizer, learning_rate, time", output_file)

    for(i in 1:3) {      ### multiple runs to increase accuracy
        mclapply(1:nrow(all_param), mc.cores = parallel_cores, function(r){
        # mclapply(1:2, mc.cores = 2, function(r){
            ro = c(all_param[r,])
            print(paste0("parameter set ",r))
            check_performance(output_file, modelType = as.character(ro$modelType), se=se, encoding_dim = ro$encoding_dim, batch_size = ro$batch_size, 
                              epochs=ro$epochs, es_patience=EARLY_STOPPING_PATIENCE, learning_rate = ro$learning_rate, seed=SEED, optimizer=OPTIMIZER, verbose = AE_VERBOSE)
        } )
    }
    
}



############################################
### creates a .csv table "hyper_par_table_full.csv" with all parameters and results 
### create new output file
output_file = file.path(RESULT_DIR, "hyperpar_table.csv")
cat(NULL,file=output_file)

run_over_param(output_file, parallel_cores = 8)   




############################################
### analyze hyperparamter runs by plot
table_hyperpar = fread(file.path(RESULT_DIR,"hyperpar_table.csv"))


### check precision recall for differnt epochs number
png(file = file.path(RESULT_DIR, "hyperpar_encoding_dim.png"), width = 600, height = 500, units = "px", res = 100)
ggplot(data=table_hyperpar, aes(x=encoding_dim_PC, y=auc)) +
  geom_smooth() +
  geom_point(shape=16, alpha=0.5) +
  labs(title ="hyperparameter optimization on inserted outliers\n[3 runs per encoding dimension]", y = "AUC precision-recall curve", x = "encoding dimension") +
  scale_color_brewer(palette="Set1") + theme(plot.title = element_text(hjust = 0.5))
dev.off()














