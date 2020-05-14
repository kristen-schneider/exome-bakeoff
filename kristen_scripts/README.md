in conda environemnt with modules:
pysam
numpy

run `python3 ./kristen_scripts/pileup_main.py -p /Shares/layer_shared/projects/chco/exome_bakeoff/mpileup_files/ -r /Shares/layer_shared/projects/chco/exome_bakeoff/Analyses/Biases/GCBias/acmg_all_59.bed -m quality`

**-p** path to downsampled pileups
**-r** path to regions direcotry (e.g. path to 59 ACMG genes)
**-m** metric which we want to calculate (e.g. quality)

## pileup_main.py


