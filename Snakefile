# DIRECTORIES
#
#

out_dir = "/Users/krsc0813/exome-bakeoff/"

rule all:
    input:
        out_dir + "hw.done"

rule hw:
    input:
        "hw.py"
    
    output:
        out_dir + 'hw.done'
    
    shell:
        "python3 {input} && touch " + out_dir + "hw.done"
