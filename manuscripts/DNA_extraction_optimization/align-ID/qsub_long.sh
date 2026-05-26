#!/bin/bash
source /etc/profile

module load nextflow
module load Python/3.11.1

config=$1

# primary nextflow execution command: alter as needed
nextflow run main.nf -c $config -profile singularity,cluster_submit,long