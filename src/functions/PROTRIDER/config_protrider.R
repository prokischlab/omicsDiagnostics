############################################
### PROTRIDER config file
############################################

### specify PROTRIDER settings
ENCODING_DIM = 25
NOISE_FACTOR = 0.5
SEED = 111 # NULL
BATCH_SIZE = 5000 # no batches during training as they increase confounders influence
EPOCHS = 500
EARLY_STOPPING_PATIENCE = 40  # early stopping patience
OPTIMIZER = "Adam"
LEARNING_RATE = 0.01
AE_VERBOSE = 1  # ae training process shown: 0 = silent, 1 = progress bar, 2 = one line per epoch

WITH_CONTROL_SAMPLES = FALSE
CONFOUNDERS_USED = c("gender", "PROTEOMICS_BATCH","BATCH_RUN", "INSTRUMENT")

### outlier limits
OUTLIER_FC = 0
OUTLIER_PVAL_ADJ = 0.05





