# this file concatenates all of teh individual acmg bed files into one bed file
# this script is responsible for creating acmg_all_59.bed

# remove the old file before appending all beds to it
rm acmg_all_59.bed
for f in /Shares/layer_shared/projects/chco/exome_bakeoff/59_genes/*.bed
do
	cat $f >> acmg_all_59.bed
done
# print out the number of genes found in acmg_all_59.bed
awk -v x=5 '{print $x}' acmg_all_59.bed | uniq | wc -l
