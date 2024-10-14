#!/bin/bash

#SBATCH --nodes=1
#SBATCH --array=1-5
#SBATCH --output=test.log
#SBATCH --job-name=run_test

# load modules
ml R/4.3

# execute script
Rscript bash_test.R $SLURM_ARRAY_TASK_ID

wait