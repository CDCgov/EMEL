#!/bin/bash
map_file=$1
fasta=$2
tax_id=$3

for i in $(grep '>' $fasta | cut -d '>' -f2 | cut -d ' ' -f1); do

    dup_check=$(cut -f1 $map_file | grep "$i")

    if [[ -z $dup_check ]]; then

        echo -e "$i\t$tax_id" >> $map_file
    
    else

        echo "Sequence ID [$i] already in mapping file [$map_file] - skipping contig"

    fi

done