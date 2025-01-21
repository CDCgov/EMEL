#!/bin/bash -l

#$ -N Extract
#$ -cwd
#$ -q all.q

## $1 == genome coordinate 1 to trim (i.e. "NW_020312409:1749-1750")
## $2 == genome coordinate 2 to trim (i.e. "NW_020312409:1874-1875")
## Extract reads that map to both coordinates
## samtools v1.13

conda activate samtools

## Remove all reads that do not contain first coordinate
for i in *_Cyclospora-mapped.bam; do
	id=$(basename $i _Cyclospora-mapped.bam)
	samtools view -h -b $i $1 > ${id}_temp.bam
done

## Reindex bam 
for i in *_temp.bam; do
	id=$(basename $i _temp.bam)
	samtools index $i
done

## Remove all reads that do not contain second coordinate
for i in *_temp.bam; do
        id=$(basename $i _temp.bam)
        samtools view -h -b $i $2 > ${id}_cut.bam
done

rm *temp*
