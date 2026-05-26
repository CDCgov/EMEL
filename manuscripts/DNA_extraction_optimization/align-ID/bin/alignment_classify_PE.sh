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
unmapped="$prefix-unmapped.sam"

echo "Parsing BAM headers"

# printing bam headers to all primary alignment output sam files
singularity exec --bind $bind_dir $SAMTOOLS_CONTAINER samtools view -H $bam > $primary_all 
singularity exec --bind $bind_dir $SAMTOOLS_CONTAINER samtools view -H $bam > $primary_unambiguous
singularity exec --bind $bind_dir $SAMTOOLS_CONTAINER samtools view -H $bam > $primary_ambiguous_single_genome
singularity exec --bind $bind_dir $SAMTOOLS_CONTAINER samtools view -H $bam > $primary_ambiguous_multi_genome
singularity exec --bind $bind_dir $SAMTOOLS_CONTAINER samtools view -H $bam > $unmapped

echo "Converting BAM to SAM with no headers"

singularity exec --bind $bind_dir $SAMTOOLS_CONTAINER samtools view $bam > "$prefix-complete.sam" #converting bam to sam for easier parsing in the loop below

complete_sam=$prefix-complete.sam


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

singularity exec --bind $bind_dir $SAMTOOLS_CONTAINER samtools view $bam | awk '$16 == "tp:A:P" || $16 == "tp:A:S" || $2 == 165 || $2 == 101 || $2 == 133 || $2 == 77 || $2 == 141 || $2 == 69' | awk '$2 < 2000' | sed 's|"|\\"|g' | awk -v OFS='\t' '{print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $16}' > $trimmed_sam

singularity exec --bind $bind_dir $SQLITE3_CONTAINER sqlite3 $sam_db \
	"
	CREATE TABLE sam_complete (read_id TEXT, flag INTEGER, ref_id TEXT, left_most_position INTEGER, mapq INTEGER, cigar TEXT, r_next TEXT, p_next TEXT, t_length TEXT, sequence TEXT, phred TEXT, nm_tag TEXT, ms_tag TEXT, as_tag TEXT, tp_tag TEXT);
	"

singularity exec --bind $bind_dir $SQLITE3_CONTAINER sqlite3 $sam_db <<EOF
.mode tabs
.import $trimmed_sam sam_complete
EOF

# paired reads which are both unmapped have their ref_id = *
# a mate that is unmapped when the other is mapped will have it's r_next = '=' and cigar = *
singularity exec --bind $bind_dir $SQLITE3_CONTAINER sqlite3 $sam_db \
	"
	UPDATE sam_complete
	SET ref_id = 'unclassified'
	WHERE ref_id = '*' or (r_next = '=' and cigar = '*')
	"

singularity exec --bind $bind_dir $SQLITE3_CONTAINER sqlite3 $sam_db \
	"
	CREATE TABLE tax_map (seq_id TEXT, tax_id TEXT);
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
	CREATE TABLE total_read_pairs (num_reads FLOAT);
	"

singularity exec --bind $bind_dir $SQLITE3_CONTAINER sqlite3 $sam_db \
	"
	INSERT INTO total_read_pairs
	SELECT COUNT(read_id) / 2
	FROM sam_complete
	WHERE tp_tag != 'tp:A:S'
	"

singularity exec --bind $bind_dir $SQLITE3_CONTAINER sqlite3 $sam_db \
	"
	CREATE TABLE unmapped_pairs (read_id TEXT, r1_flag INTEGER, r2_flag INTEGER, ref_id TEXT, mapq INTEGER);
	"

singularity exec --bind $bind_dir $SQLITE3_CONTAINER sqlite3 $sam_db \
	"
	INSERT INTO unmapped_pairs
	SELECT read_id, MIN(flag), MAX(flag), ref_id, mapq
	FROM sam_complete
	WHERE ref_id = 'unclassified' and (flag = 77 or flag = 141)
	GROUP BY read_id
	"

singularity exec --bind $bind_dir $SQLITE3_CONTAINER sqlite3 $sam_db \
	"
	CREATE TABLE r1_reads (read_id TEXT, r1_flag INTEGER, r1_ref_id TEXT, r1_tax_id TEXT, r1_mapq INTEGER, r1_p_next INTEGER, r1_as_tag INTEGER);
	"

singularity exec --bind $bind_dir $SQLITE3_CONTAINER sqlite3 $sam_db \
	"
	INSERT INTO r1_reads
	SELECT read_id, MIN(flag), ref_id, tax_map.tax_id, mapq, p_next, as_tag
	FROM sam_complete
	LEFT JOIN tax_map ON (sam_complete.ref_id = tax_map.seq_id)
	WHERE tp_tag = 'tp:A:P' or (r_next = '=' and cigar = '*')
	GROUP BY read_id
	"

singularity exec --bind $bind_dir $SQLITE3_CONTAINER sqlite3 $sam_db \
	"
	CREATE TABLE r2_reads (read_id TEXT, r2_flag INTEGER, r2_ref_id TEXT, r2_tax_id TEXT, r2_mapq INTEGER, r2_p_next INTEGER, r2_as_tag INTEGER);
	"

singularity exec --bind $bind_dir $SQLITE3_CONTAINER sqlite3 $sam_db \
	"
	INSERT INTO r2_reads
	SELECT read_id, MAX(flag), ref_id, tax_map.tax_id, mapq, p_next, as_tag
	FROM sam_complete
	LEFT JOIN tax_map ON (sam_complete.ref_id = tax_map.seq_id)
	WHERE tp_tag = 'tp:A:P' or (r_next = '=' and cigar = '*')
	GROUP BY read_id
	"

singularity exec --bind $bind_dir $SQLITE3_CONTAINER sqlite3 $sam_db \
	"
	CREATE TABLE mapped_pairs (read_id TEXT, r1_flag INTEGER, r2_flag INTEGER, r1_ref_id TEXT, r2_ref_id TEXT, r1_tax_id TEXT, r2_tax_id TEXT, r1_mapq INTEGER, r2_mapq INTEGER, r1_p_next INTEGER, r2_p_next INTEGER, r1_as_tag INTEGER, r2_as_tag INTEGER);
	"

singularity exec --bind $bind_dir $SQLITE3_CONTAINER sqlite3 $sam_db \
	"
	INSERT INTO mapped_pairs
	SELECT r1.read_id, r1.r1_flag, r2.r2_flag, r1.r1_ref_id, r2.r2_ref_id, r1.r1_tax_id, r2.r2_tax_id, r1.r1_mapq, r2.r2_mapq, r1.r1_p_next, r2.r2_p_next, r1.r1_as_tag, r2.r2_as_tag
	FROM r1_reads as r1
	JOIN r2_reads as r2 ON (r1.read_id = r2.read_id)
	GROUP BY r1.read_id
	"

# mapped_pairs table may contain read pairs where one mate is unmapped while the other is mapped - here we are converting their reference ID and tax ID to unclassified
singularity exec --bind $bind_dir $SQLITE3_CONTAINER sqlite3 $sam_db \
	"
	UPDATE mapped_pairs
	SET r1_tax_id = 'unclassified'
	WHERE r1_ref_id = 'unclassified' and r1_tax_id IS NULL;

	UPDATE mapped_pairs
	SET r2_tax_id = 'unclassified'
	WHERE r2_ref_id = 'unclassified' and r2_tax_id IS NULL;
	"

singularity exec --bind $bind_dir $SQLITE3_CONTAINER sqlite3 $sam_db \
	"
	CREATE TABLE top_secondary_mappings (read_id TEXT, sec_tax_id TEXT, sec_p_next INTEGER, sec_as_tag INTEGER);
	"

singularity exec --bind $bind_dir $SQLITE3_CONTAINER sqlite3 $sam_db \
	"
	INSERT INTO top_secondary_mappings
	SELECT read_id, tax_map.tax_id as tax_id, p_next, MAX(as_tag) as as_tag
	FROM sam_complete
	LEFT JOIN tax_map ON (sam_complete.ref_id = tax_map.seq_id)
	WHERE tp_tag = 'tp:A:S'
	GROUP BY read_id, p_next, tax_id
	"

singularity exec --bind $bind_dir $SQLITE3_CONTAINER sqlite3 $sam_db \
	"
	CREATE TABLE top_secondary_mappings_taxa_count (read_id TEXT, taxa_count INTEGER);
	"

singularity exec --bind $bind_dir $SQLITE3_CONTAINER sqlite3 $sam_db \
	"
	INSERT INTO top_secondary_mappings_taxa_count
	SELECT read_id, COUNT(DISTINCT sec_tax_id) as taxa_count
	FROM top_secondary_mappings
	GROUP BY read_id
	"

singularity exec --bind $bind_dir $SQLITE3_CONTAINER sqlite3 $sam_db \
	"
	CREATE TABLE secondary_mapping_aggregate (read_id TEXT, sec_tax_id TEXT, sec_p_next INTEGER, sec_as_tag INTEGER, sec_taxa_count INTEGER);
	"

singularity exec --bind $bind_dir $SQLITE3_CONTAINER sqlite3 $sam_db \
	"
	INSERT INTO secondary_mapping_aggregate
	SELECT top_secondary_mappings.read_id, top_secondary_mappings.sec_tax_id, top_secondary_mappings.sec_p_next, top_secondary_mappings.sec_as_tag, top_secondary_mappings_taxa_count.taxa_count
	FROM top_secondary_mappings
	LEFT JOIN top_secondary_mappings_taxa_count ON (top_secondary_mappings.read_id = top_secondary_mappings_taxa_count.read_id)
	"

singularity exec --bind $bind_dir $SQLITE3_CONTAINER sqlite3 $sam_db \
	"
	CREATE TABLE unambiguous_pairs (read_id TEXT, r1_flag INTEGER, r2_flag INTEGER, r1_ref_id TEXT, r2_ref_id TEXT, r1_tax_id TEXT, r2_tax_id TEXT, r1_mapq INTEGER, r2_mapq INTEGER, r1_p_next INTEGER, r2_p_next INTEGER, r1_as_tag INTEGER, r2_as_tag INTEGER);
	"

singularity exec --bind $bind_dir $SQLITE3_CONTAINER sqlite3 $sam_db \
	"
	INSERT INTO unambiguous_pairs
	SELECT *
	FROM mapped_pairs
	WHERE r1_tax_id = r2_tax_id and r1_mapq > 0 and r2_mapq > 0
	"

singularity exec --bind $bind_dir $SQLITE3_CONTAINER sqlite3 $sam_db \
	"
	CREATE TABLE ambiguous_single_genome_pairs (read_id TEXT, r1_tax_id TEXT, r2_tax_id TEXT, r1_mapq INTEGER, r2_mapq INTEGER, sec_tax_id TEXT, sec_taxa_count INTEGER);
	"

singularity exec --bind $bind_dir $SQLITE3_CONTAINER sqlite3 $sam_db \
	"
	INSERT INTO ambiguous_single_genome_pairs
	SELECT mp.read_id, mp.r1_tax_id, mp.r2_tax_id, mp.r1_mapq, mp.r2_mapq, sma.sec_tax_id, sma.sec_taxa_count
	FROM mapped_pairs as mp
	JOIN secondary_mapping_aggregate as sma ON mp.read_id = sma.read_id AND mp.r1_tax_id = sma.sec_tax_id
	WHERE mp.r1_tax_id = mp.r2_tax_id and (mp.r1_mapq = 0 or mp.r2_mapq = 0) and sma.sec_taxa_count = 1
	GROUP BY mp.read_id
	
	UNION ALL

	SELECT mp.read_id, mp.r1_tax_id, mp.r2_tax_id, mp.r1_mapq, mp.r2_mapq, sma.sec_tax_id, sma.sec_taxa_count
	FROM mapped_pairs as mp
	JOIN secondary_mapping_aggregate as sma ON mp.read_id = sma.read_id AND mp.r1_tax_id = sma.sec_tax_id
	WHERE mp.r1_tax_id != 'unclassified' and mp.r2_tax_id = 'unclassified' and mp.r1_mapq = 0 and sma.sec_taxa_count = 1
	GROUP BY mp.read_id

	UNION ALL

	SELECT mp.read_id, mp.r1_tax_id, mp.r2_tax_id, mp.r1_mapq, mp.r2_mapq, sma.sec_tax_id, sma.sec_taxa_count
	FROM mapped_pairs as mp
	JOIN secondary_mapping_aggregate as sma ON mp.read_id = sma.read_id AND mp.r2_tax_id = sma.sec_tax_id
	WHERE mp.r1_tax_id = 'unclassified' and mp.r2_tax_id != 'unclassified' and mp.r2_mapq = 0 and sma.sec_taxa_count = 1
	GROUP BY mp.read_id

	UNION ALL

	SELECT mp.read_id, mp.r1_tax_id, mp.r2_tax_id, mp.r1_mapq, mp.r2_mapq, sma.sec_tax_id, sma.sec_taxa_count
	FROM mapped_pairs as mp
	JOIN secondary_mapping_aggregate as sma ON mp.read_id = sma.read_id AND mp.r2_tax_id = sma.sec_tax_id
	WHERE mp.r1_tax_id != 'unclassified' and mp.r2_tax_id = 'unclassified' and mp.r1_mapq > 0
	GROUP BY mp.read_id
		
	UNION ALL

	SELECT mp.read_id, mp.r1_tax_id, mp.r2_tax_id, mp.r1_mapq, mp.r2_mapq, sma.sec_tax_id, sma.sec_taxa_count
	FROM mapped_pairs as mp
	JOIN secondary_mapping_aggregate as sma ON mp.read_id = sma.read_id AND mp.r2_tax_id = sma.sec_tax_id
	WHERE mp.r1_tax_id = 'unclassified' and mp.r2_tax_id != 'unclassified' and mp.r2_mapq > 0
	GROUP BY mp.read_id
	"

singularity exec --bind $bind_dir $SQLITE3_CONTAINER sqlite3 $sam_db \
	"
	CREATE TABLE ambiguous_multi_genome_pairs (read_id TEXT, r1_tax_id TEXT, r2_tax_id TEXT, r1_mapq INTEGER, r2_mapq INTEGER, r1_as_tag INTEGER, r2_as_tag INTEGER);
	"

singularity exec --bind $bind_dir $SQLITE3_CONTAINER sqlite3 $sam_db \
	"
	INSERT INTO ambiguous_multi_genome_pairs
	SELECT mp.read_id, mp.r1_tax_id, mp.r2_tax_id, mp.r1_mapq, mp.r2_mapq, mp.r1_as_tag, mp.r2_as_tag
	FROM mapped_pairs as mp
	WHERE mp.read_id NOT IN (
		SELECT read_id 
		FROM unambiguous_pairs

		UNION ALL

		SELECT read_id
		FROM ambiguous_single_genome_pairs
	)
	"

singularity exec --bind $bind_dir $SQLITE3_CONTAINER sqlite3 $sam_db \
	"
	CREATE TABLE unambiguous_hits_by_taxid (prim_tax_id TEXT, total_hits);
	"

singularity exec --bind $bind_dir $SQLITE3_CONTAINER sqlite3 $sam_db \
	"
    INSERT INTO unambiguous_hits_by_taxid
	SELECT r1_tax_id, COUNT(DISTINCT read_id)
	FROM unambiguous_pairs
	GROUP BY r1_tax_id;
	"

singularity exec --bind $bind_dir $SQLITE3_CONTAINER sqlite3 $sam_db \
	"
	CREATE TABLE ambiguous_single_genome_hits_by_taxid (prim_tax_id TEXT, total_hits);
	"

# for single genome ambiguous pairs we need to make sure the primary tax ID is set to the tax ID of the mate that is NOT 'unclassified' when applicable
singularity exec --bind $bind_dir $SQLITE3_CONTAINER sqlite3 $sam_db \
	"
	INSERT INTO ambiguous_single_genome_hits_by_taxid
	SELECT q1.prim_tax_id, COUNT(q1.prim_tax_id)
	FROM (
	SELECT CASE
			WHEN r1_tax_id = 'unclassified' THEN
				r2_tax_id
			WHEN r2_tax_id = 'unclassified' THEN
				r1_tax_id
			ELSE
				r1_tax_id
			END prim_tax_id
	FROM ambiguous_single_genome_pairs
	) as q1
	GROUP BY q1.prim_tax_id
	"

singularity exec --bind $bind_dir $SQLITE3_CONTAINER sqlite3 $sam_db \
	"
	CREATE TABLE ambiguous_multi_genome_hits_by_taxid (prim_tax_id TEXT, total_hits);
	"

singularity exec --bind $bind_dir $SQLITE3_CONTAINER sqlite3 $sam_db \
	"
	INSERT INTO ambiguous_multi_genome_hits_by_taxid
	SELECT q1.prim_tax_id, COUNT(q1.prim_tax_id)
	FROM (
	SELECT CASE
			WHEN r1_tax_id = 'unclassified' THEN
				r2_tax_id
			WHEN r2_tax_id = 'unclassified' THEN
				r1_tax_id
			WHEN r1_as_tag > r2_as_tag THEN
				r1_tax_id
			WHEN r2_as_tag > r1_as_tag THEN
				r2_tax_id			
			ELSE
				r1_tax_id
			END prim_tax_id
	FROM ambiguous_multi_genome_pairs
	) as q1
	GROUP BY q1.prim_tax_id
	"

echo -e "species\tnum_reads\tpercent_reads\tambiguity" > "$summary_file"

singularity exec --bind $bind_dir $SQLITE3_CONTAINER sqlite3 $sam_db \
	"
    SELECT unmapped_pairs.ref_id, COUNT(unmapped_pairs.read_id), ROUND(COUNT(unmapped_pairs.read_id) / total_read_pairs.num_reads * 100,6) as percentage
    FROM unmapped_pairs
	LEFT JOIN total_read_pairs
	GROUP BY ref_id
	" | sed 's/|/\t/g' | awk -v OFS='\t' '{print $0, "none"}' >> "$summary_file"

singularity exec --bind $bind_dir $SQLITE3_CONTAINER sqlite3 $sam_db \
	"
    SELECT unambiguous_hits_by_taxid.prim_tax_id, unambiguous_hits_by_taxid.total_hits, ROUND(unambiguous_hits_by_taxid.total_hits / total_read_pairs.num_reads * 100,6) as percentage
    FROM unambiguous_hits_by_taxid
	LEFT JOIN total_read_pairs
	GROUP BY prim_tax_id
	" | sed 's/|/\t/g' | awk -v OFS='\t' '{print $0, "none"}' >> "$summary_file"

singularity exec --bind $bind_dir $SQLITE3_CONTAINER sqlite3 $sam_db \
	"
    SELECT ambiguous_single_genome_hits_by_taxid.prim_tax_id, ambiguous_single_genome_hits_by_taxid.total_hits, ROUND(ambiguous_single_genome_hits_by_taxid.total_hits / total_read_pairs.num_reads * 100,6) as percentage
    FROM ambiguous_single_genome_hits_by_taxid
	LEFT JOIN total_read_pairs
	GROUP BY prim_tax_id
	" | sed 's/|/\t/g' | awk -v OFS='\t' '{print $0, "single_genome"}' >> "$summary_file"

singularity exec --bind $bind_dir $SQLITE3_CONTAINER sqlite3 $sam_db \
	"
    SELECT ambiguous_multi_genome_hits_by_taxid.prim_tax_id, ambiguous_multi_genome_hits_by_taxid.total_hits, ROUND(ambiguous_multi_genome_hits_by_taxid.total_hits / total_read_pairs.num_reads * 100,6) as percentage
    FROM ambiguous_multi_genome_hits_by_taxid
	LEFT JOIN total_read_pairs
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
	
	elif [[ "$tax_id" =~ [a-zA-Z] ]]; then # if the tax_id var contains letters then it is most likely a nonstandard reference (or it wasn't included in the taxonomy map) and we don't want to try and look this up in the names.dmp file

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
	WHERE tp_tag = 'tp:A:P' and read_id IN (
		SELECT read_id
		FROM mapped_pairs
	)
	" | sed 's/!##########!/\t/g' >> $primary_all

# grab all nonambiguous primary alignments
singularity exec --bind $bind_dir $SQLITE3_CONTAINER sqlite3 $sam_db ".separator '!##########!'" \
	"
    SELECT *
    FROM sam_complete
	WHERE tp_tag = 'tp:A:P' and read_id IN (
		SELECT read_id
		FROM unambiguous_pairs
	)
	" | sed 's/!##########!/\t/g' >> $primary_unambiguous

# grab all primary ambiguous single genome hits
singularity exec --bind $bind_dir $SQLITE3_CONTAINER sqlite3 $sam_db ".separator '!##########!'" \
	"
    SELECT *
    FROM sam_complete
	WHERE tp_tag = 'tp:A:P' and read_id IN (
		SELECT read_id
		FROM ambiguous_single_genome_pairs
	)
	" | sed 's/!##########!/\t/g' >> $primary_ambiguous_single_genome

# grab all primary ambiguous multi-genome hits
singularity exec --bind $bind_dir $SQLITE3_CONTAINER sqlite3 $sam_db ".separator '!##########!'" \
	"
    SELECT *
    FROM sam_complete
	WHERE tp_tag = 'tp:A:P' and read_id IN (
		SELECT read_id
		FROM ambiguous_multi_genome_pairs
	) 
	" | sed 's/!##########!/\t/g' >> $primary_ambiguous_multi_genome

# grab all unmapped reads
singularity exec --bind $bind_dir $SQLITE3_CONTAINER sqlite3 $sam_db ".separator '!##########!'" \
	"
    SELECT *
    FROM sam_complete
	WHERE read_id IN (
		SELECT read_id
		FROM unmapped_pairs
	) 
	" | sed 's/!##########!/\t/g' >> $unmapped


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

singularity exec --bind $bind_dir $SAMTOOLS_CONTAINER samtools view -@ 16 -b $unmapped > $prefix-unmapped.bam
singularity exec --bind $bind_dir $SAMTOOLS_CONTAINER samtools sort $prefix-unmapped-sorted.bam > $prefix-unmapped-sorted.bam
rm -f $prefix-unmapped.bam

# cleaning up temporary files
rm -f $prefix-chunk_*
rm -f $prefix-part-* 
rm -f sam_chunk*
rm -f *.sam
rm -f $prefix-ambig.*