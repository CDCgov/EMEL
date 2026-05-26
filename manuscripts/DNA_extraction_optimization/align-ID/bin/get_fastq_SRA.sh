#!/bin/bash

for i in $@; do

    # retrieve fastq files from SRA
    wget -O $i.fastq.gz "https://www.ncbi.nlm.nih.gov/Traces/sra-reads-be/fastq?acc=$i&clipped=1"

done