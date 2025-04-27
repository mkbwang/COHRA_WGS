#!/bin/bash
# The interpreter used to execute the script

#“#SBATCH” directives that convey submission options:

#SBATCH --job-name=plaque_batch_correction
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --nodes=1
#SBATCH --cpus-per-task=12
#SBATCH --mem-per-cpu=6g
#SBATCH --time=02:00:00
#SBATCH --account=bfoxman0
#SBATCH --partition=standard
#SBATCH --output=logs/plaque_batch_correction.log

module load Rtidyverse/4.2.0

Rscript batch_correction_plaque.R

