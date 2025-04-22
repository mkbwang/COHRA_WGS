#!/bin/bash

#set input directory
in=concat_reads/

 
# get count of files in this directory
NUMFILES=$(ls -1 ${in}*_cat.fq.gz | wc -l)

# subtract 1 as we have to use zero-based indexing (first element is 0)
ZBNUMFILES=$(($NUMFILES - 1))
#ZBNUMFILES=196

# submit array of jobs to SLURM
if [ $ZBNUMFILES -ge 0 ]; then
  sbatch --array=0-${ZBNUMFILES} array_humann.pbs
else
  echo "No jobs to submit, since no input files in this directory."
fi
