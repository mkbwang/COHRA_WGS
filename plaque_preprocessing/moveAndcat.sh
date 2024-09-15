#! /bin/bash
#moveAndcat.sh
#F Blostein
#University of Michigan
# modified by Mukai Wang
#Purpose: Move trimmed, cleaned and decontaminated files from individual sample folders to another folder and concatonate them
#for use with humann2
#Usage: bash moveAndcat.sh

##########################
#Set Script Env#
##########################

#Variables 
my_dir=paired_reads/
# my_dir=/nfs/turbo/sph-bfoxman/mcbyrd/MetaWGS/Cleaning/Post_Decontam/Post_Host_Removal/
#new_dir=/scratch/bfoxman_root/bfoxman/blostein/saliva_meta/biobakery_saliva/cat_files
new_dir=concat_reads/

cd ${my_dir}
for basename in $(ls *.fq.gz | cut -f1 -d_ | sort -u)
do
    cat "$basename"*R1*.fq.gz "$basename"*R2*.fq.gz > "../${new_dir}/${basename}_cat.fq.gz"
done

end_files=$(ls ${new_dir} | wc -l)
start_files=$(ls ${my_dir}*R1* | wc -l)

echo "started with $start_files and ended with $end_files"
