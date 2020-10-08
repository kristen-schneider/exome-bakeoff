#!/usr/bin/env bash
#
#SBATCH -p short
#SBATCH --job-name=intersect_regions_ss
#SBATCH --ntasks=1
#SBATCH --time=4:00:00
#SBATCH --mem-per-cpu=10G
#SBATCH --output=/Shares/layer_shared/projects/chco/exome_bakeoff/regions/splice_sites/ss.out
#SBATCH --error=/Shares/layer_shared/projects/chco/exome_bakeoff/regions/splice_sites/ss.err

# PURPOSE: intersect all splice sites with all bed files from a given region

regions='/Shares/layer_shared/projects/chco/exome_bakeoff/regions/oncogenes/'
splice_sites='/Shares/layer_shared/projects/chco/exome_bakeoff/regions/splice_sites/ss.bed'
outdir='/Shares/layer_shared/projects/chco/exome_bakeoff/regions/splice_sites/'

echo $regions

for f_bed in `ls $regions`; do
    if [[ $f_bed == *.bed ]]; then
        echo $f_bed
        
        # prepare the name of the output file by removing the pattern '.*' greedily
        sample=${f_bed%%.*}
        out_f=${outdir}${sample}_ss.bed

        # Check if the output file already exists, if so, skip the conversion
        if test -e $out_f; then
            echo "$out_f already exists, so skipping processing it"

        else
            echo "submitting a job"
            qsub intersect_worker.sh $regions$f_bed $splice_sites $out_f
        fi
    fi
    
done

#exons_bed=/Shares/layer_shared/projects/chco/exome_bakeoff/regions/Homo_sapiens.GRCh37.82.exons.bed
#SCAP_vcf=/Shares/layer_shared/projects/chco/exome_bakeoff/regions/intersect_bed_vcf/scap_COMBINED_v1.0.vcf

#exons_bed=/Shares/layer_shared/projects/chco/exome_bakeoff/regions/intersect_bed_vcf/test.bed
#SCAP_vcf=/Shares/layer_shared/projects/chco/exome_bakeoff/regions/intersect_bed_vcf/test.vcf

#bedtools intersect  -a $exons_bed -b $SCAP_vcf > output_intersect.bed
