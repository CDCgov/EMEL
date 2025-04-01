#!/bin/bash -l

#$ -N Align-and-sort
#$ -cwd
#$ -q all.q

#### $1 == indexed reference genome
#### Map to reference using strict parameters
#### bowtie2 v2.3.5.1
#### samtools v1.13 

module load bowtie2
conda activate samtools

cd Test_dataset_fastq/

for i in *.fastq; do
	id=$(basename $i .fastq)
	bowtie2 -p 4 -x $1 -U $i --score-min L,0,-0.2 | \
	samtools view -bu | samtools sort -@6 -l 6 - \
	> ${id}.bam
done

sh ../bam-statistics.sh

mkdir BAM
mv *.bam BAM

