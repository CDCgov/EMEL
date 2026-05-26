#!/bin/bash
R_CONTAINER="https://depot.galaxyproject.org/singularity/r-basejump%3A0.18.1--r44hdfd78af_0"

bind_dir="/"$(echo $(pwd) | cut -d '/' -f2)

prefix=$1
alignment_report=$2
rscript=$3
data_type=$4

input_tsv="$prefix-rscript-input.csv"

singularity exec --bind $bind_dir $R_CONTAINER Rscript $rscript $alignment_report "$(pwd)" $prefix $data_type

mv plot.png $prefix-taxonomic-distribution-plot_mqc.png