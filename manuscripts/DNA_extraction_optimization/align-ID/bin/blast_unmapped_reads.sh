#!/bin/bash
SAMTOOLS_CONTAINER="https://depot.galaxyproject.org/singularity/samtools%3A1.21--h50ea8bc_0"
BLAST_CONTAINER="https://depot.galaxyproject.org/singularity/blast%3A2.15.0--pl5321h6f7f691_1"

bind_dir="/"$(echo $(pwd) | cut -d '/' -f2)

prefix=$1
bam=$2
blastdb=$3
blast_evalue=$4
blast_perc_identity=$5
blast_target_seqs=$6

singularity exec --bind $bind_dir $SAMTOOLS_CONTAINER samtools fasta $bam > $prefix-unmapped.fasta

random_fasta_ids=$(cat $prefix-unmapped.fasta | grep '>' | shuf -n 50)

IFS=$'\n'

for i in $random_fasta_ids; do

    grep -A 1 "$i" $prefix-unmapped.fasta >> $prefix-unmapped-subsample.fasta

done

BLASTDB=/scicomp/groups-pure/WDPB/EMEL/Projects/Long_Read_Analysis_RUSHER/data/blast/taxonomy
blast_args="6 qseqid qstart qend sseqid sstart send pident length mismatch evalue bitscore qcovhsp qcovs sacc stitle staxids"

singularity exec --bind $bind_dir $BLAST_CONTAINER blastn \
    -db "$blastdb" \
    -query $prefix-unmapped-subsample.fasta \
    -outfmt "$blast_args" \
    -num_threads 16 \
    -evalue $blast_evalue \
    -perc_identity $blast_perc_identity \
    -max_target_seqs $blast_target_seqs \
    -qcov_hsp_perc 90 \
    -out $prefix-unmapped.txt

mod_blast_args="${blast_args/6 /}"
mod_blast_args=$(echo "$mod_blast_args" | sed 's/ /\t/g')
echo -e "$mod_blast_args" > $prefix-unmapped-blast-hits_sorted.txt

cat $prefix-unmapped.txt | sort -k12,12 -nr >> $prefix-unmapped-blast-hits_sorted.txt

rm -f $prefix-unmapped.txt