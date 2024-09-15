#!/bin/sh

# make sure you have activated the biobakery environment

source /home/wangmk/.bashrc
micromamba activate biobakery

DIR=metaphlan

# transform the counts from SGB to GTDB
allfiles=(ls $DIR/*_cat.txt)
for file in "${allfiles[@]}"; do
    base="${file%.*}"
    sgb_to_gtdb_profile.py -i ${DIR}/${file} -o ${DIR}/${base}_gtdb.txt
done

# merge the counts of individual files
merge_metaphlan_tables.py ${DIR}/*_cat.txt > ${DIR}/joined_taxonomic_profile.tsv
merge_metaphlan_tables.py ${DIR}/*_cat_gtdb.txt > ${DIR}/joined_taxonomic_profile_gtdb.tsv

# identify all the unique species identified among all the samples
humann_reduce_table --input ${DIR}/joined_taxonomic_profile.tsv --output ${DIR}/max_taxonomic_profile.tsv --function max --sort-by level

grep -e "s__" ${DIR}/max_taxonomic_profile.tsv | grep -v "t__"  > ${DIR}/merged_taxonomic_table_species.tsv

cut -d "|" -f 6-7 ${DIR}/merged_taxonomic_table_species.tsv > ${DIR}/merged_taxonomic_table_species_concise.tsv

sed -i '1i# mpa_vJun23_CHOCOPhlAnSGB_202307' ${DIR}/merged_taxonomic_table_species_concise.tsv

sed -i 's/$/\t /' ${DIR}/merged_taxonomic_table_species_concise.tsv
