#!/bin/bash

# sbatch -N 1 -c 12 --mem 20G --time 04:00:00 --error job.%J.err --output job.%J.out -J hyperopt -w ouga04 hyperopt_ouga04.sh

#sbatch -N 1 -c 15 --mem 20G --error job.%J.err --output job.%J.out -J hyperopt -w ouga04 hyperopt_run.sh

sbatch -N 1 -c 8 --mem 10G --error job.%J.err --output job.%J.out -J hyperopt hyperopt_run.sh




#srun -c 12 --mem=64G --pty Rscript hyperparameter_optimization_complete.R
                  


