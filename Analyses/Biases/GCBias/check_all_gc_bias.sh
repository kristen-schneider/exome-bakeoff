cd /scratch/Shares/layer/projects/chco/bake_off_temp_files/bams_filtered_to_acmg_59/
for f in *.bam
do
	name="$(echo $f | sed -e 's/.bqsr.bam//g')"
	echo $name
	dir="/Shares/layer_shared/projects/chco/exome_bakeoff/Analyses/Biases/GCBias/results/"
	dir="$dir$name"
	#echo $dir
	mkdir $dir
	input="/scratch/Shares/layer/projects/chco/bake_off_temp_files/bams_filtered_to_acmg_59/"
	input="$input$f"
	output="$dir/gc_bias_metrics.txt"
	chart="$dir/gc_bias_metrics.pdf"
	summary="$dir/summary_metrics.txt"
	bash /Shares/layer_shared/projects/chco/exome_bakeoff/Analyses/Biases/GCBias/gc_bias_metrics.sh $input $output $chart $summary
done
