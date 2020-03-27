cd results
for d in *
do
	file="$d/gc_bias_metrics.pdf"
	echo $file
	newfile="/Users/mibr6115/GCBiasResults/$d"
	extension="_gc_bias_metrics.pdf"
	cp $file "$newfile$extension"
done
