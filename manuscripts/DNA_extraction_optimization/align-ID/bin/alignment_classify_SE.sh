#!/bin/bash
SAMTOOLS_CONTAINER="https://depot.galaxyproject.org/singularity/samtools%3A1.21--h50ea8bc_0"
SQLITE3_CONTAINER="https://depot.galaxyproject.org/singularity/sqlite%3A3"

bind_dir="/"$(echo $(pwd) | cut -d '/' -f2)

prefix=$1
bam=$2
seqid2taxid=$3
filter_alignment_by_id=$4
my_tax_ids=$5
include_children=$6
nodes=$7
taxa_names=$8
project_dir=$9
nodes_sqlite=${10}
non_standard_reference=${11}
mapq=${12}
single_end=${13}

# Determining child nodes based on parent tax ID provided
if [[ "$filter_alignment_by_id" == "true" && "$include_children" == "true" ]]; then

    echo "Collecting all children tax ids from SQL database"

    >"$prefix-parent_and_child_ids.txt"

    for i in $(cat $my_tax_ids); do

        data=$(singularity exec --bind $bind_dir $SQLITE3_CONTAINER sqlite3 $nodes_sqlite "SELECT * FROM TAX_IDS WHERE parent_id = "$i"")

        parent_id="$i"        
        echo $parent_id >> $prefix-parent_and_child_ids.txt

        if [[ -n $data ]]; then

            child_ids=$(echo $data | cut -d '|' -f3 | sed 's/,/ /g' )

            for i in $child_ids; do

                echo $i >> $prefix-parent_and_child_ids.txt
            
            done
        
        fi

    done

    tax_ids_i_want="$prefix-parent_and_child_ids.txt"

else

    tax_ids_i_want="placeholder"

fi

# if a sam file was used as input, we need to convert it to bam format first
if [[ $bam == *".sam" ]]; then

    singularity exec --bind $bind_dir $SAMTOOLS_CONTAINER samtools sort -@ 4 -o "$prefix".bam $bam
    bam="$prefix".bam

fi

# defining all alignment output files
summary_file="$prefix-taxonomic-summary.tsv"
primary_all="$prefix-primary_all.sam"
primary_unambiguous="$prefix-primary_unambiguous.sam"
primary_ambiguous_single_genome="$prefix-primary_ambiguous_single_genome.sam"
primary_ambiguous_multi_genome="$prefix-primary_ambiguous_multi_genome.sam"

echo "Parsing BAM headers"

# printing bam headers to all primary alignment output sam files
singularity exec --bind $bind_dir $SAMTOOLS_CONTAINER samtools view -H $bam > $primary_all 
singularity exec --bind $bind_dir $SAMTOOLS_CONTAINER samtools view -H $bam > $primary_unambiguous
singularity exec --bind $bind_dir $SAMTOOLS_CONTAINER samtools view -H $bam > $primary_ambiguous_single_genome
singularity exec --bind $bind_dir $SAMTOOLS_CONTAINER samtools view -H $bam > $primary_ambiguous_multi_genome

echo "Converting BAM to SAM with no headers"

singularity exec --bind $bind_dir $SAMTOOLS_CONTAINER samtools view $bam > "$prefix-complete.sam" #converting bam to sam for easier parsing in the loop below

complete_sam=$prefix-complete.sam

# grabbing all unmapped reads
singularity exec --bind $bind_dir $SAMTOOLS_CONTAINER samtools view -f 4 -h $bam > "$prefix-unmapped.sam"
singularity exec --bind $bind_dir $SAMTOOLS_CONTAINER samtools view -@ 8 -S -b "$prefix-unmapped.sam" > "$prefix-unmapped.bam"

if [[ "$filter_alignment_by_id" == "true" ]]; then

    tax_filtered_sam="$prefix-tax-filtered.sam"
    singularity exec --bind $bind_dir $SAMTOOLS_CONTAINER samtools view -H $bam > "$tax_filtered_sam" #printing bam headers to output sam file

fi


#############################################################################################
###                                 SAM Parsing Section                                  ####
#############################################################################################

sam_db="$prefix-sam.db"
trimmed_sam="$prefix-trimmed.sam"

rm -f $sam_db


# trimming SAM file for entries we want, the 'sed' line will escape any rogue double quotes ["] in the phred score - if we don't do this, SQLite won't be able to import the SAM file
singularity exec --bind $bind_dir $SAMTOOLS_CONTAINER samtools view $bam | awk '$16 == "tp:A:P" || $16 == "tp:A:S" || $2 == 4' | awk '$2 != 2048 && $2 != 2064' | sed 's|"|\\"|g' | awk -v OFS='\t' '{print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $16}' > $trimmed_sam


singularity exec --bind $bind_dir $SQLITE3_CONTAINER sqlite3 $sam_db \
	"
	CREATE TABLE sam_complete (read_id TEXT, flag INTEGER, ref_id TEXT, left_most_position INTEGER, mapq INTEGER, cigar TEXT, r_next TEXT, p_next INTEGER, t_length TEXT, sequence TEXT, phred BLOB, nm_tag TEXT, ms_tag TEXT, as_tag TEXT, tp_tag TEXT);
	"

singularity exec --bind $bind_dir $SQLITE3_CONTAINER sqlite3 $sam_db <<EOF
.mode tabs
.import $trimmed_sam sam_complete
EOF


singularity exec --bind $bind_dir $SQLITE3_CONTAINER sqlite3 $sam_db \
	"
	UPDATE sam_complete
	SET ref_id = 'unclassified'
	WHERE ref_id = '*' and flag = 4
	"

singularity exec --bind $bind_dir $SQLITE3_CONTAINER sqlite3 $sam_db \
	"
	CREATE TABLE tax_map (seq_id TEXT, tax_id INTEGER);
	"

singularity exec --bind $bind_dir $SQLITE3_CONTAINER sqlite3 $sam_db <<EOF
.mode tabs
.import $seqid2taxid tax_map
EOF

# manually adding reference sequence IDs to taxonomy map if they had not already been added
singularity exec --bind $bind_dir $SQLITE3_CONTAINER sqlite3 $sam_db \
	"
	INSERT INTO tax_map
	SELECT ref_id, ref_id
	FROM sam_complete
	WHERE ref_id NOT IN (
		SELECT seq_id
		FROM tax_map
	)
	"

# removing duplicate rows from seq id to tax id map file
singularity exec --bind $bind_dir $SQLITE3_CONTAINER sqlite3 $sam_db \
	"
	DELETE FROM tax_map
	WHERE rowid NOT IN (
		SELECT MIN(rowid)
		FROM tax_map
		GROUP BY seq_id	
	)
	"

singularity exec --bind $bind_dir $SQLITE3_CONTAINER sqlite3 $sam_db \
	"
	CREATE TABLE total_reads (num_reads FLOAT);
	"

singularity exec --bind $bind_dir $SQLITE3_CONTAINER sqlite3 $sam_db \
	"
	INSERT INTO total_reads
	SELECT COUNT(read_id)
	FROM sam_complete
	WHERE flag = 4 or flag = 16 or flag = 0;
	"

singularity exec --bind $bind_dir $SQLITE3_CONTAINER sqlite3 $sam_db \
	"
	CREATE TABLE unmapped (read_id TEXT, flag INTEGER, ref_id TEXT, mapq INTEGER);
	"

singularity exec --bind $bind_dir $SQLITE3_CONTAINER sqlite3 $sam_db \
	"
	INSERT INTO unmapped
	SELECT read_id, flag, ref_id, mapq
	FROM sam_complete
	WHERE ref_id = 'unclassified'
	"

singularity exec --bind $bind_dir $SQLITE3_CONTAINER sqlite3 $sam_db \
	"
	CREATE TABLE primary_multimap (read_id TEXT, flag INTEGER, ref_id TEXT, mapq INTEGER, as_tag TEXT, tax_id INTEGER, FOREIGN KEY (ref_id) REFERENCES tax_map (seq_id));
	"

singularity exec --bind $bind_dir $SQLITE3_CONTAINER sqlite3 $sam_db \
	"
	INSERT INTO primary_multimap
	SELECT sam_complete.read_id, sam_complete.flag, sam_complete.ref_id, sam_complete.mapq, sam_complete.as_tag, tax_map.tax_id
	FROM sam_complete
	JOIN tax_map ON (sam_complete.ref_id = tax_map.seq_id)
	WHERE tp_tag = 'tp:A:P' and mapq = 0
	"

singularity exec --bind $bind_dir $SQLITE3_CONTAINER sqlite3 $sam_db \
	"
	CREATE TABLE secondary_mappings (sec_read_id TEXT, sec_flag INTEGER, sec_ref_id TEXT, sec_mapq INTEGER, sec_as_tag TEXT, sec_tax_id INTEGER, FOREIGN KEY (sec_read_id) REFERENCES primary_multimap (read_id), FOREIGN KEY (sec_ref_id) REFERENCES tax_map (seq_id));
	"

singularity exec --bind $bind_dir $SQLITE3_CONTAINER sqlite3 $sam_db \
	"
	INSERT INTO secondary_mappings
	SELECT sam_complete.read_id, sam_complete.flag, sam_complete.ref_id, sam_complete.mapq, sam_complete.as_tag, tax_map.tax_id
	FROM sam_complete
	JOIN tax_map ON (sam_complete.ref_id = tax_map.seq_id)
	WHERE tp_tag = 'tp:A:S' and read_id IN (SELECT read_id FROM primary_multimap)
	"

singularity exec --bind $bind_dir $SQLITE3_CONTAINER sqlite3 $sam_db \
	"
	CREATE TABLE primary_secondary_aggregate (read_id TEXT, sec_tax_id INTEGER, sec_tax_id_count INTEGER, prim_tax_id INTEGER);
	"

singularity exec --bind $bind_dir $SQLITE3_CONTAINER sqlite3 $sam_db \
	"
	INSERT INTO primary_secondary_aggregate
	SELECT secondary_mappings.sec_read_id, secondary_mappings.sec_tax_id, COUNT(secondary_mappings.sec_tax_id), primary_multimap.tax_id
	FROM secondary_mappings
	JOIN primary_multimap ON (secondary_mappings.sec_read_id = primary_multimap.read_id) 
	WHERE secondary_mappings.sec_as_tag = primary_multimap.as_tag or cast(substr(secondary_mappings.sec_as_tag, 6) AS INTEGER) > cast(substr(primary_multimap.as_tag, 6) AS INTEGER)
	GROUP BY secondary_mappings.sec_read_id, secondary_mappings.sec_tax_id
	"

singularity exec --bind $bind_dir $SQLITE3_CONTAINER sqlite3 $sam_db \
	"
	CREATE TABLE primary_unambiguous (read_id TEXT, flag INTEGER, ref_id TEXT, mapq INTEGER, as_tag TEXT, tax_id INTEGER, FOREIGN KEY (ref_id) REFERENCES tax_map (seq_id));
	"

# had to change the condition for WHERE clause - a MAPQ of 0 doesn't guarantee a secondary alignment with a >= alignment score (AS) as the primary alignment meaning that is read may not be accounted for in the ambiguous (single or multi genome) queries; so in order to grab all unambiguous primary alignments we will use the tp_tag and make sure the read_id does not appear in the primary_secondary_aggregate table
singularity exec --bind $bind_dir $SQLITE3_CONTAINER sqlite3 $sam_db \
	"
	INSERT INTO primary_unambiguous
	SELECT sam_complete.read_id, sam_complete.flag, sam_complete.ref_id, sam_complete.mapq, sam_complete.as_tag, tax_map.tax_id
	FROM sam_complete
	JOIN tax_map ON (sam_complete.ref_id = tax_map.seq_id)
	WHERE sam_complete.tp_tag = 'tp:A:P' AND sam_complete.read_id NOT IN (
		SELECT read_id
		FROM primary_secondary_aggregate
	)
	"

singularity exec --bind $bind_dir $SQLITE3_CONTAINER sqlite3 $sam_db \
	"
	CREATE TABLE single_genome_hits_by_read (read_id TEXT, sec_tax_id INTEGER, sec_tax_id_count INTEGER, prim_tax_id INTEGER);
	"	

singularity exec --bind $bind_dir $SQLITE3_CONTAINER sqlite3 $sam_db \
	"
	INSERT INTO single_genome_hits_by_read
	SELECT read_id, sec_tax_id, sec_tax_id_count, prim_tax_id
	FROM primary_secondary_aggregate
	WHERE sec_tax_id = prim_tax_id 
	GROUP BY read_id, prim_tax_id
	HAVING COUNT(read_id) = 1
	"

singularity exec --bind $bind_dir $SQLITE3_CONTAINER sqlite3 $sam_db \
	"
	CREATE TABLE single_genome_hits_by_taxid (prim_tax_id INTEGER, total_hits);
	"

singularity exec --bind $bind_dir $SQLITE3_CONTAINER sqlite3 $sam_db \
	"
	INSERT INTO single_genome_hits_by_taxid
	SELECT prim_tax_id, COUNT(prim_tax_id)
	FROM single_genome_hits_by_read
	GROUP BY prim_tax_id
	" 

singularity exec --bind $bind_dir $SQLITE3_CONTAINER sqlite3 $sam_db \
	"
	CREATE TABLE multi_genome_hits_by_read (read_id TEXT, sec_tax_id INTEGER, sec_tax_id_count INTEGER, prim_tax_id INTEGER);
	"

singularity exec --bind $bind_dir $SQLITE3_CONTAINER sqlite3 $sam_db \
	"
	INSERT INTO multi_genome_hits_by_read
	SELECT read_id, sec_tax_id, sec_tax_id_count, prim_tax_id
	FROM primary_secondary_aggregate
	WHERE read_id NOT IN (
		SELECT read_id
		FROM single_genome_hits_by_read
	) 
	" 

singularity exec --bind $bind_dir $SQLITE3_CONTAINER sqlite3 $sam_db \
	"
	CREATE TABLE multi_genome_hits_by_taxid (prim_tax_id INTEGER, total_hits);
	"

singularity exec --bind $bind_dir $SQLITE3_CONTAINER sqlite3 $sam_db \
	"
	INSERT INTO multi_genome_hits_by_taxid
	SELECT prim_tax_id, COUNT(DISTINCT read_id)
	FROM multi_genome_hits_by_read
	GROUP BY prim_tax_id
	"

singularity exec --bind $bind_dir $SQLITE3_CONTAINER sqlite3 $sam_db \
	"
	CREATE TABLE unambiguous_hits_by_taxid (prim_tax_id INTEGER, total_hits);
	"

singularity exec --bind $bind_dir $SQLITE3_CONTAINER sqlite3 $sam_db \
	"
    INSERT INTO unambiguous_hits_by_taxid
	SELECT tax_id, COUNT(DISTINCT read_id)
	FROM primary_unambiguous
	GROUP BY tax_id;
	"

echo -e "species\tnum_reads\tpercent_reads\tambiguity" > "$summary_file"

singularity exec --bind $bind_dir $SQLITE3_CONTAINER sqlite3 $sam_db \
	"
    SELECT unmapped.ref_id, COUNT(unmapped.read_id), ROUND(COUNT(unmapped.read_id) / total_reads.num_reads * 100,6) as percentage
    FROM unmapped
	LEFT JOIN total_reads
	GROUP BY ref_id
	" | sed 's/|/\t/g' | awk -v OFS='\t' '{print $0, "none"}' >> "$summary_file"

singularity exec --bind $bind_dir $SQLITE3_CONTAINER sqlite3 $sam_db \
	"
    SELECT unambiguous_hits_by_taxid.prim_tax_id, unambiguous_hits_by_taxid.total_hits, ROUND(unambiguous_hits_by_taxid.total_hits / total_reads.num_reads * 100,6) as percentage
    FROM unambiguous_hits_by_taxid
	LEFT JOIN total_reads
	GROUP BY prim_tax_id
	" | sed 's/|/\t/g' | awk -v OFS='\t' '{print $0, "none"}' >> "$summary_file"

singularity exec --bind $bind_dir $SQLITE3_CONTAINER sqlite3 $sam_db \
	"
    SELECT single_genome_hits_by_taxid.prim_tax_id, single_genome_hits_by_taxid.total_hits, ROUND(single_genome_hits_by_taxid.total_hits / total_reads.num_reads * 100,6) as percentage
    FROM single_genome_hits_by_taxid
	LEFT JOIN total_reads
	GROUP BY prim_tax_id
	" | sed 's/|/\t/g' | awk -v OFS='\t' '{print $0, "single_genome"}' >> "$summary_file"

singularity exec --bind $bind_dir $SQLITE3_CONTAINER sqlite3 $sam_db \
	"
    SELECT multi_genome_hits_by_taxid.prim_tax_id, multi_genome_hits_by_taxid.total_hits, ROUND(multi_genome_hits_by_taxid.total_hits / total_reads.num_reads * 100,6) as percentage
    FROM multi_genome_hits_by_taxid
	LEFT JOIN total_reads
	GROUP BY prim_tax_id
	" | sed 's/|/\t/g' | awk -v OFS='\t' '{print $0, "multi_genome"}' >> "$summary_file"

count=0

while IFS= read -r line; do

	#skipping header line
	if [ $count -lt 1 ]; then

		count=$((count+1))
		continue
	
	fi

	tax_id=$(echo "$line" | cut -f1)

	if [[ $tax_id == "unclassified" ]]; then

		continue
	
	elif [[ "$tax_id" =~ [a-zA-Z] ]]; then # if the tax_id var contains letters then it is most likely a nonstandard reference

		continue

	fi

	# this line is grabbing the scientific name from the names.dmp taxonomy file
	species_name=$(grep "^$tax_id.[|].*scientific name" $taxa_names | cut -d "|" -f 2 | tr -d '\t')

	sed -i "s|^\<$tax_id\>|$species_name|g" "$summary_file"

done < "$summary_file"


# grab all primary alignments
singularity exec --bind $bind_dir $SQLITE3_CONTAINER sqlite3 $sam_db ".separator '!##########!'" \
	"
    SELECT *
    FROM sam_complete
	WHERE tp_tag = 'tp:A:P'
	" | sed 's/!##########!/\t/g' >> $primary_all

# grab all nonambiguous primary alignments
singularity exec --bind $bind_dir $SQLITE3_CONTAINER sqlite3 $sam_db ".separator '!##########!'" \
	"
    SELECT *
    FROM sam_complete
	WHERE tp_tag = 'tp:A:P' and read_id IN (
		SELECT read_id
		FROM primary_unambiguous
	)
	" | sed 's/!##########!/\t/g' >> $primary_unambiguous

# grab all primary ambiguous single genome hits
singularity exec --bind $bind_dir $SQLITE3_CONTAINER sqlite3 $sam_db ".separator '!##########!'" \
	"
    SELECT *
    FROM sam_complete
	WHERE tp_tag = 'tp:A:P' and mapq = 0 and read_id IN (
		SELECT read_id
		FROM single_genome_hits_by_read
	)
	" | sed 's/!##########!/\t/g' >> $primary_ambiguous_single_genome

# grab all primary ambiguous multi-genome hits
singularity exec --bind $bind_dir $SQLITE3_CONTAINER sqlite3 $sam_db ".separator '!##########!'" \
	"
    SELECT *
    FROM sam_complete
	WHERE tp_tag = 'tp:A:P' and mapq = 0 and read_id NOT IN (
		SELECT read_id
		FROM single_genome_hits_by_read
	) 
	" | sed 's/!##########!/\t/g' >> $primary_ambiguous_multi_genome


#############################################################################################
#############################################################################################
#############################################################################################


if [[ "$filter_alignment_by_id" == "true" ]]; then

    #converting sam file into bam file
    singularity exec --bind $bind_dir $SAMTOOLS_CONTAINER samtools view -@ 16 -S -b $tax_filtered_sam > $prefix-classified-plus-filtered.bam

    #sorting BAM file
    singularity exec --bind $bind_dir $SAMTOOLS_CONTAINER samtools sort $prefix-classified-plus-filtered.bam > $prefix-classified-plus-filtered-sorted.bam

fi


singularity exec --bind $bind_dir $SAMTOOLS_CONTAINER samtools view -@ 16 -b $primary_all > $prefix-primary-all.bam
singularity exec --bind $bind_dir $SAMTOOLS_CONTAINER samtools sort $prefix-primary-all.bam > $prefix-primary-all-sorted.bam
rm -f $prefix-primary-all.bam

singularity exec --bind $bind_dir $SAMTOOLS_CONTAINER samtools view -@ 16 -b $primary_unambiguous > $prefix-primary-unambiguous.bam
singularity exec --bind $bind_dir $SAMTOOLS_CONTAINER samtools sort $prefix-primary-unambiguous.bam > $prefix-primary-unambiguous-sorted.bam
rm -f $prefix-primary-unambiguous.bam

singularity exec --bind $bind_dir $SAMTOOLS_CONTAINER samtools view -@ 16 -b $primary_ambiguous_single_genome > $prefix-primary_ambiguous_single_genome.bam
singularity exec --bind $bind_dir $SAMTOOLS_CONTAINER samtools sort $prefix-primary_ambiguous_single_genome.bam > $prefix-primary_ambiguous_single_genome-sorted.bam
rm -f $prefix-primary_ambiguous_single_genome.bam

singularity exec --bind $bind_dir $SAMTOOLS_CONTAINER samtools view -@ 16 -b $primary_ambiguous_multi_genome > $prefix-primary_ambiguous_multi_genome.bam
singularity exec --bind $bind_dir $SAMTOOLS_CONTAINER samtools sort $prefix-primary_ambiguous_multi_genome.bam > $prefix-primary_ambiguous_multi_genome-sorted.bam
rm -f $prefix-primary_ambiguous_multi_genome.bam

# cleaning up temporary files
rm -f $prefix-chunk_*
rm -f $prefix-part-* 
rm -f sam_chunk*
rm -f *.sam
rm -f $prefix-ambig.*