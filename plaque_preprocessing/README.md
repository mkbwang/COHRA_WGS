This folder contains multiple `.sh` and `.pbs` scripts. The majority of the taxonomy annotation and gene function annotation were done with [metaphlan4](https://huttenhower.sph.harvard.edu/metaphlan/) and [humann3](https://huttenhower.sph.harvard.edu/humann/) in the biobakery pipeline. The uses of them are:

1. `count_seq.sh`: count the number fo reads in each fastq file.
2. `array_decontam.pbs`and `metagenomeDecontaminateReads.sh` and `decontam.sh`: Remove the reads that map to spike-in species genomes (Truepera radiovictrix, Allobacillus halotolerans, Imtechella halotolerans)
3. `move_cat.pbs` and `moveAndcat.sh`: concatenate the forward reads and backward reads of the same fastq file.
4. `array_metaphlan.pbs` and `metaphlan.sh`: Run metaphlan to get the taxa counts of each sample.
5. `merge_metaphlan.sh`: Merge the relative abundance of each sample into one file.
6. `single_humann.pbs`: Run humann for one sample. This step includes setting up the reference nucleotide database according to species identified by metaphlan.
7. `array_humann.pbs` and `humann.sh`: Run humann for all the rest of the samples.
8. `humann_merge_regroup.pbs`: Merge the pathway/gene abundance of each single file into one file.
