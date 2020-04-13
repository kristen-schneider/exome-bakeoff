txt_dir='/Shares/layer_shared/projects/chco/exome_bakeoff/pileup-tbx/txt/'
output_dir='/Shares/layer_shared/projects/chco/exome_bakeoff/pileup-tbx/gz/'
bed_gz_extension='.bed.gz'

pip install pysam

for txt in `ls ${txt_dir}`
do
    # find base name. convert to bed. bgzip.
    IFS='.' # space is set as delimiter
    read -ra AD <<< "$txt" # str is read into an array as tokens separated by IFS
    cat "$txt_dir$txt" | awk '{OFS="\t"; print $1,$2,$2+1,$3,$4,$5,$6,$7;}' | bgzip -c >  "$output_dir$AD$bed_gz_extension" 
done
