#!/bin/bash
# The interpreter used to execute the script

#“#SBATCH” directives that convey submission options:

#SBATCH --job-name=saliva_subset_genes
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=10g
#SBATCH --time=00:10:00
#SBATCH --account=bfoxman0
#SBATCH --partition=standard
#SBATCH --output=logs/saliva_subset_genes.log

reffile="/nfs/turbo/sph-bfoxman/People/wangmk/COHRA_WGS/saliva_preprocessing/saliva_valid_unirefs.txt"
inputfolder="/nfs/turbo/sph-bfoxman/People/wangmk/COHRA_WGS/saliva_preprocessing/humann_output"
outputfolder="/nfs/turbo/sph-bfoxman/People/wangmk/COHRA_WGS/saliva_preprocessing/humann_output/subset_genes"


inputfiles=($(find $inputfolder -maxdepth 1 -type f -name "*genefamilies.tsv"))

for fname in ${inputfiles[@]}
do 
  file=($(basename ${fname}))
  grep -F -f $reffile $fname > ${outputfolder}/${file}
done

