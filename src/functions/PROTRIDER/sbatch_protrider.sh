#!/bin/bash
#SBATCH --time=20:00:00
#SBATCH --nodes=1
#SBATCH --nodelist=ouga04
#SBATCH --cpus-per-task=20
#SBATCH --job-name=protrider
#SBATCH --output=py.%j.out
#SBATCH --error=py.%j.err
#SBATCH --mem=15G


### 1. clone https://gitlab.cmm.in.tum.de/gagneurlab/py_outrider
### 2. specify PYTHONPATH below
### 3. switch to an environment with required python packages installed: conda activate loipfi_tf2
### 4. run this script in command line: sbatch sbatch_protrider.sh
### 5. edit output list to fit column names to previous ones (link to RSCRIPT)


export PYTHONPATH="/data/nasif12/home_if12/smirnovd/multiOMICs_integration/py_outrider/py_outrider/"

filepath='/s/project/mitoMultiOmics/proteome_analysis/'
cpu=20
covariates="PROTEOMICS_BATCH BATCH_RUN gender INSTRUMENT"


python /data/nasif12/home_if12/smirnovd/multiOMICs_integration/py_outrider/py_outrider/main.py --file_meas $filepath"/processed_data/protrider_intensities_filtered_input.csv" --profile protrider --output "/s/project/mitoMultiOmics/multiOMICs_integration/processed_data/protrider/py_output" --num_cpus=$cpu --max_iter 15 --file_sa $filepath"/processed_data/protrider_proteomics_annotation.csv" --covariates $covariates --output_plots True --output_list True --seed 5





