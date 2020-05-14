# General Workflow
in conda environemnt with modules:
pysam
numpy

run the following line of code, first <br>
`python3 ./kristen_scripts/pileup_main.py`<br>
  `    -p /Shares/layer_shared/projects/chco/exome_bakeoff/mpileup_files/`<br>
  `    -r /Shares/layer_shared/projects/chco/exome_bakeoff/Analyses/Biases/GCBias/acmg_all_59.bed`<br>
  `    -m quality`<br>

**-p** path to downsampled pileups<br>
**-r** path to regions direcotry (e.g. path to 59 ACMG genes)<br>
**-m** metric which we want to calculate (e.g. quality)<br>

# Specific file functions
### pileup_main.py
iterates through the directory of pileup files. for each pileup file, it will call the script single_pileup.py.

### sigle_pileup.py



