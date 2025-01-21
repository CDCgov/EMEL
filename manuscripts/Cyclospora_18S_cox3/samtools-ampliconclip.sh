#!/bin/bash -l

#$ -N Extract
#$ -cwd
#$ -q all.q

## $1 == bedfile
## Clips overhangs of mapped reads in bam files based on bedfile coordinates
## samtools v1.13 

module load conda
conda activate samtools

for i in *cut.bam; do
	id=$(basename $i _cut.bam)
	samtools ampliconclip -o ${id}_clipped.bam \
	-b $1  $i --both-ends --hard-clip
done

