dir='fastq_downsample'
for x in `ls ${dir}/*_R1_001.fq`
do
    bn=`basename $x`
    sample=${bn%%_R1_001.fq}
    echo -e "$sample \t $x \t ${dir}/${sample}_R2_001.fq"
    #ls $dir/${sample}_R2_001.fq
done
