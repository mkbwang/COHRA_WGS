#!/bin/bash


#SBATCH --job-name=stats_plaque_fastq
#SBATCH--mail-type=BEGIN,FAIL,END
#SBATCH --nodes=1
#SBATCH --cpus-per-task=6
#SBATCH --mem=15g
#SBATCH --time=6:00:00
#SBATCH --account=bfoxman0
#SBATCH --partition=standard
#SBATCH --output=logs/stats_fastq_plaque.log

source /home/wangmk/.bashrc
micromamba activate biobakery

#input_folder="/nfs/turbo/sph-bfoxman/Sequences/PostHR"
#seqkit -j 4 stats ${input_folder}/*.fq.gz > stats_with_spikein.tsv

output_folder="./fastq_files"
seqkit -j 6 stats ${output_folder}/*.fq.gz > stats_without_spikein.tsv

