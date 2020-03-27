# $1 should be the bam file you want to filter
# $2 should be the output location for the filtered bam
echo "Working on:"
echo $1
echo "----"
samtools view -h -b -L /Shares/layer_shared/projects/chco/exome_bakeoff/Analyses/Biases/GCBias/acmg_all_59.bed $1 > $2
