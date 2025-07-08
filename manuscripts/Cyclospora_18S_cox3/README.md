# [Cyclospora_18S_COX3](https://github.com/CDCgov/WDPB_EMEL/tree/main/manuscripts/Cyclospora_18S_cox3)
### Associated Publication
#### [Hofstetter, J., Arfken, A., Kahler, A., Qvarnstrom, Y., Rodrigues, C., & Mattioli, M. (2024). Evaluation of coccidia DNA in irrigation pond water and wastewater sludge associated with Cyclospora cayetanensis 18S rRNA gene qPCR detections. Microbiology Spectrum, e00906-24.](https://journals.asm.org/doi/full/10.1128/spectrum.00906-24)

## Project Description
Selected scripts and codes used in the dectection of *Cyclospora cayetanensis* from 18S and COX3 amplicon sequencing. Test data sets included in this repository are from a subset of samples pulled from COX3 sequencing.


## Bash Scripts
`Blast-NCBI-nt.sh` use ncbi-blast+ with fasta files to blast against a selected database with a 90% ID & coverage cutoff and return results with tab-delimited format and the following parameters: queryID accession taxID title length percentID query coverage mismatch gapopen query start query end subject start subject end evalue bitscore.   
          
`samtools-extract-coordinates.sh` extract bam alignments based on genome coordinates from indexed bam files. Requires two coordinates. Test with Test_dataset_GENESLICE.
    
Example Usage:
```
sh samtools-extract-coordinates.sh "NW_020312409:1749-1750" "NW_020312409:1874-1875"
```
`samtools-ampliconclip.sh` trim overhangs from extracted bam alignments. Requires bed file with designated trimming locations. Test with Test_dataset_GENESLICE following coordinate extraction.   
     
Example Usage:
```
sh samtools-ampliconclip.sh 1749-1875-bed
```
## Python codes

`Parse_Blast_Dict_dupcounts.py` parse blast outputs and duplicate sequence count data from seqkit rmdup (https://github.com/shenwei356/seqkit) and return table with combined top counts based on max bitscores for all samples. A randomly selected hit will be selected if max bitscores are tied among several hits. Note: tied hits may include all hits with the same taxID, does not distinguish. Requires a blast output file and count data for each sample. Use with Test_dataset_BLAST.
    
An optional taxa `-target` may also be selected to retrieve all max bitscore hits based on taxa ID. Note: for target option, tied hits are only considered tied maxscore hits between different taxID (i.e tied hits among same taxID will considered unique and not tied). **** *Warning: inconsistent counts may occur when using the string option if taxa have different names or are missing labels for the same taxID!* ****
  
General Parameters:
```
-d          directory containing files
-s1         file suffix for blast output results
-s2         file suffix for duplicate count files
-o          prefix for output files (optional)
-target     character string or taxID (optional)
```
   
Example Usage:
```
python Parse_Blast_Dict_dupcounts.py -d . -s1 _blast-hits-nt.txt -s2 _duplicated.detail.txt
```
Example Usage with Targets:
```
python Parse_Blast_Dict_dupcounts.py -d . -s1 _blast-hits-nt.txt -s2 _duplicated.detail.txt -o Cyclo_ -target 88456
 ```
 ```
python Parse_Blast_Dict_dupcounts.py -d . -s1 _blast-hits-nt.txt -s2 _duplicated.detail.txt -o Eimeria_ -target Eimeria
 ```
    
`FastqList.py` create a list of representative sequence IDs for each sample based on a community percent threshold. Use with Test_dataset_FASTQLIST.

Parameters:
```
-d          directory containing files
-s          file suffix of duplicate count files
-p          threshold percentage (given as a decimal)
```
    
Example Usage:
```
python FastqList.py -d . -s _duplicated.detail.txt -p 0.25
```
## Figures
`18S_Dendrogram.md` rmarkdown to create dendrogram plot of most abundant taxa found in pond and sludge samples based on 18S gene sliced reads mapped to *C. cayetanensis* 18S gene region (655-807bp).

## Metadata
Org: Division of Foodborne Waterborne and Environmental Diseases, WDPB

Version: 1.0.0 

Status: Maintained

Keywords: pathogen detection, next generation sequencing analysis, bioinformatics, genomics, metagenomics, parasitic disease, *Cyclospora cayetanensis*

Contact Email: ncezid_shareit@cdc.gov

