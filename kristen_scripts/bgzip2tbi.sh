bgzip_dir='/Shares/layer_shared/projects/chco/exome_bakeoff/pileup-tbx/'
output_dir='/Shares/layer_shared/projects/chco/exome_bakeoff/pileup-tbx/'

pip install pysam

for gz in `ls ${bgzip_dir}`
do
    # tabix
    tabix "$bgzip_dir$gz" 
done

