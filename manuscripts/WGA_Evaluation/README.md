# Salmonella Detection Limit (WGA)


## Project Description
Selected scripts and codes used for determining the detection limit and sensitivity of  *Salmonella* .spp in mixed microbially communities using whole meta-/genome amplification (WMA/WGA). 
    
Test data sets included in this repository are from a subset of samples pulled from the ATCC Environmental mock community dataset. Test depth and fastq files are subsampled at depth of 250,000 reads each (merged paired-end).
    
## Bash Scripts
`Align-and-stats.sh` map fastq reads to reference using strict mapping threshold and calls `bam-statistics.sh` to return depth, coverage, and mapping statistics. Set to run in Test_dataset_fastq directory.
    
Example Usage:
```
sh Align-and-stats.sh ../Salmonella_REF/Salmonella
```
    
## Python codes
`Parse_Blast_custom_WGA.py` parse BLAST output from FastQ-screen multihit reads using customized ATCC Environmental database and ncbi-blast+. The following format and delimiters are required in the outputs to be parsed: 
    
-outfmt "6 qseqid sacc staxid stitle length pident qcovs mismatch gapopen qstart qend sstart send evalue bitscore" \  

Parameters:
```
-d      directory
-s      file suffix for blast output results
-mock   mock community (Environmental, Zymo, Isolate)
```
    
Use with Test_dataset_multihit.

Example Usage:
```
python Parse_Blast_custom_WGA.py -d Test_dataset_multihit -s _blast-hits.txt -mock Environmental
```
   
`CV.py` return table of calculated coefficient of variation in genome depths based on samtools depth profile outputs. Use with Test_dataset_depth.

Parameters:
```
-d      directory
-s      file suffix for depth files
-r      gene/contig/plasmid/region on genome
```
Eample Usage:
```
python CV.py -d Test_dataset_depth -s _depth.txt -r NC_003197.2
```
    
## Jupyter Notebooks

`Depth_plots.ipynb` create visual depth plots from samtools depth profile outputs. Use with Test_dataset_depth.

## Mansucript status
In preparation
