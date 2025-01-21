#!/bin/bash -l

#$ -N blast
#$ -cwd
#q all.q
#$ -pe smp 12

#### Blast query template script for nt scicomp nt database

module load ncbi-blast+

for i in *.fasta; do
	id=$(basename $i .fasta)
	blastn -db /scicomp/reference/ncbi-blast/current/nt \
	-perc_identity 90 \
	-qcov_hsp_perc 90 \
	-max_hsps 1 \
	-max_target_seqs 5 \
	-outfmt "6 qseqid sacc staxid stitle length pident qcovs mismatch gapopen qstart qend sstart send evalue bitscore" \
	-num_threads $NSLOTS \
	-query $i \
	-out ${id}_blast-hits-nt.txt
done

