#!/bin/bash
#SBATCH -p general
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=2g
#SBATCH -t 2-5:00:00

module add matlab
matlab -nodesktop -nosplash -singleCompThread -r dyn_code -logfile mycode.out
