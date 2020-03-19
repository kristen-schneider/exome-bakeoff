# for every bam file in the following dir
# use samtools to select only the information relating to acmg 59
# saves the filtered files in /scratch/Shares/layer/projects/chco/bake_off_temp_files/bams_filtered_to_acmg_59/
cd /Shares/layer_shared/projects/chco/exome_bakeoff/bam_files
for f in *.bam
do
	fullfile="/Shares/layer_shared/projects/chco/exome_bakeoff/bam_files/"
	fullfile="$fullfile$f"
	outfile="/scratch/Shares/layer/projects/chco/bake_off_temp_files/bams_filtered_to_acmg_59/"
	outfile="$outfile$f"
	/Shares/layer_shared/projects/chco/exome_bakeoff/Analyses/Biases/GCBias/filter_single_bam.sh $fullfile $outfile 
done
wait
