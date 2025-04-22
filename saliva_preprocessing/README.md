This folder contains multiple `.sh` and `.pbs` scripts. The majority of the taxonomy annotation and gene function annotation were done with [metaphlan4](https://huttenhower.sph.harvard.edu/metaphlan/) and [humann3](https://huttenhower.sph.harvard.edu/humann/) in the biobakery pipeline. The uses of them are:


1. `array_decontam.pbs`and `metagenomeDecontaminateReads.sh` : Remove the reads that map to spike-in species genomes (Truepera radiovictrix, Allobacillus halotolerans, Imtechella halotolerans). For saliva this is redundant since there were no spike-ins.
2. `move_cat.pbs` and `moveAndcat.sh`: concatenate the forward reads and backward reads of the same fastq file.
3. `array_metaphlan.pbs` and `metaphlan.sh`: Run metaphlan to get the taxa counts of each sample. `merge_metaphlan.sh`: Merge the relative abundance of each sample into one file.
4. `single_humann.pbs`: Run humann for one sample and get uniref90 protein abundances. This step includes setting up the reference nucleotide database according to species identified by metaphlan. `array_humann.pbs` and `humann.sh`: Run humann for all the rest of the samples.
5. Some uniref90s are no longer considered valid proteins in the uniprot database. I look up all the unique uniref90 proteins in the database (`uniref_split.py`, `uniref_filter_1.py`, `uniref_filter_2.py`)
6. Use `subset_genes.sh` to remove invalid genes. Then use `humann_merge_regroup.pbs` to merge the pathway/gene abundance of each single file into one file. I also merge the uniref90s into KEGG orthologs.

