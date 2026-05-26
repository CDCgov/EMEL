# Read Mapping Ambiguity with Examples

## Summary
Align-ID generates classification reports that group reads together by species and their overall level of ambiguity. The ambiguity level can fall into one of three categories: `none`, `single_genome`, or `multi-genome`. Read aligners, such as minimap2, will attempt to align each read to a location in the reference genome with the best alignment score. Sometimes a read will map equally well to more the one location in the reference genome(s), each mapping location recieving and equal alignment score. These types of alignments are considered ambiguous and the read aligner itself will randomly assign one ambiguous mapping location to be the primary alignment and any other equally scoring alignment to be secondary. In such cases in can be useful to know how many reads were ambiguously aligned and if they aligned equally as well to genomic regions within the same reference genome/species or if equal scoring alignments were generated across multiple reference geneomes/species. 

## Long and Single-End Read Ambiguity Levels
Single-end and long read ambiguity is relatively straight forward. See below for a definition of each ambiguity level:
- `none`: the read is either unmapped or aligns to only one location/species
- `single_genome`: the read aligns equally well to several location within the same reference genome/species
- `multi-genome`: the read aligns eqaully well to several locations across multiple reference genomes/species

## Paired-End Ambiguity Levels
Paired-End alignments are a bit trickier, so we'll go through each ambiguity level with examples:

### Non-ambiguous Alignments
Non-ambiguous alignments will recieve the `none` ambiguity classifier:

![Non-ambiguous alignment scenarios](images/ambiguity_none.jpg)

### Single-Genome Ambiguity
Alignments that align to multiple locations within the same refence genome or across reference genomes that belong to the same species alignments will recieve the `single-genome` ambiguity classifier:

![Single-Genome ambiguous alignment scenario 1 and 2](images/single_genome_ambiguity_1.jpg)
![Single-Genome ambiguous alignment scenario 3](images/single_genome_ambiguity_2.jpg)

### Multi-Genome Ambiguity
Alignments that align to multiple locations across more than one species of reference genomes will recieve the `multi-genome` ambiguity classifier:

![Multi-Genome ambiguous alignment scenario 1](images/multi_genome_ambiguity_1.jpg)
![Multi-Genome ambiguous alignment scenario 2](images/multi_genome_ambiguity_2.jpg)
![Multi-Genome ambiguous alignment scenario 3](images/multi_genome_ambiguity_3.jpg)
![Multi-Genome ambiguous alignment scenario 4](images/multi_genome_ambiguity_4.jpg)
