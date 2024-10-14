#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=30
#SBATCH --array=1-100
#SBATCH --output=PSA.log
#SBATCH --job-name=run_psa
#SBATCH --mem=60G
#SBATCH --time=47:59:00
#SBATCH --mail-type=ALL

# --ntasks-per-node will be used to specify the number
# of cores to use on the machine.

# working directory 
cd /home/users/mbkauf/OlderKidneyTransplantSim/

# load modules
ml R/4.3

# execute script
Rscript R/04_PSA.R $SLURM_ARRAY_TASK_ID

wait
