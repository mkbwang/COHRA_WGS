#!/bin/sh

# make sure you have activated the biobakery environment

source /home/wangmk/.bashrc
micromamba activate biobakery

DIR=metaphlan

merge_metaphlan_tables.py ${DIR}/*_cat.txt > ${DIR}/joined_taxonomic_profile.tsv

humann_reduce_table --input ${DIR}/joined_taxonomic_profile.tsv --output ${DIR}/max_taxonomic_profile.tsv --function max --sort-by level

grep -e "s__" ${DIR}/max_taxonomic_profile.tsv | grep -v "t__"  > ${DIR}/merged_taxonomic_table_species.tsv

cut -d "|" -f 6-7 ${DIR}/merged_taxonomic_table_species.tsv > ${DIR}/merged_taxonomic_table_species_concise.tsv

sed -i '1i# mpa_vOct22_CHOCOPhlAnSGB_202212' ${DIR}/merged_taxonomic_table_species_concise.tsv

sed -i 's/$/\t /' ${DIR}/merged_taxonomic_table_species_concise.tsv
