#!/bin/bash -l

#$ -N depth
#$ -cwd
#$ -q all.q

module load conda
conda activate samtools

for i in *.bam; do
	id=$(basename $i .bam)
	samtools depth -a $i > ${id}_depth.txt
done

mkdir samtools_depth
mv *_depth.txt samtools_depth

for i in *.bam; do
	id=$(basename $i .bam)
	samtools coverage $i > ${id}_coverage.txt
done

mkdir samtools_coverage
mv *_coverage.txt samtools_coverage/

for i in *.bam; do 
	id=$(basename $i .bam)
	samtools flagstat $i > ${id}_mapping.txt
done

mkdir samtools_mapping
mv *_mapping.txt samtools_mapping
