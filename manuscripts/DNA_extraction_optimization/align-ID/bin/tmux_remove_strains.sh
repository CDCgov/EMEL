#!/bin/bash

# This script divides a taxonomy nodes.dmp file into n chunks and processes each chunk in a separate tmux session.

# Get the directory where the script is located
SCRIPTDIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

# Define script params
while getopts ":f:n:s:c:h" opt; do
  case $opt in
    f) nodes_file="$OPTARG" ;;
    n) num_groups="$OPTARG" ;;
    c) seq_tax_id_map_copy="$OPTARG" ;;
    h) help=true ;;
    \?) echo "Invalid option -$OPTARG" >&2
        exit 1
    ;;
    :) echo "Option -$OPTARG requires an argument." >&2
       exit 1
    ;;
  esac
done

if [[ -z $nodes_file || -z $num_groups || -z $seq_tax_id_map_copy ]]; then
  echo -e "Usage: $0\n\n-f <PATH>\tpath to nodes.dmp taxonomy file\n-n <INT>\tspecify the number of batches/tmux sessions to split up the nodes.dmp file into - the lines in nodes.dmp will be evenly distributed amongst batches\n-c <PATH> separate copy of seq id to tax id file to be overwritten (all strain tax ids will be converted to the parent species tax id). The seqid2taxid.map file can usually be found in a kraken2 database. You should use either the standard kraken2 database or large, curaited database of similar size for best results. NOTE:\n-h\t\tprint this message"
  exit 1
fi

# Clearing temp directory
rm -rf $SCRIPTDIR/../temp/*

# creating output directories
mkdir $SCRIPTDIR/../temp/logs
mkdir $SCRIPTDIR/../temp/segments
mkdir $SCRIPTDIR/../temp/ids

LOGFILE="$SCRIPTDIR/../temp/logs/log.txt"


# Clearing log file
>$LOGFILE

SHUFFLED_LIST="$SCRIPTDIR/../temp/shuffled_list"

# Get the total number of lines in the input file
total_lines=$(wc -l $nodes_file | cut -d ' ' -f 1)
line_count=0

# Calculate the number of lines per chunk (divided into 10 chunks)
div=$((total_lines/$num_groups))
start=1
fin=$div

>$SHUFFLED_LIST
sort -R $nodes_file > $SHUFFLED_LIST

# Loop through 10 chunks
for i in $(seq 1 $num_groups); do

    # Define the input file for this chunk
    nodes_segement="$SCRIPTDIR/../temp/segments/lines_$start-$fin"
    >$nodes_segement

    # define output file for this chunk
    output_file="$SCRIPTDIR/../temp/ids/lines_$start-$fin"'_'"species-strain-ids"
    >$output_file

    # Extract lines from the input file into the output file
    sed -n ''"$start"','"$fin"'p' $SHUFFLED_LIST > $nodes_segement

    # Define a range of lines for this chunk
    lines="$start-$fin"

    # Create a tmux session for processing this chunk
    session_name="strain_replace_$start-$fin"
    echo -e "\nNew tmux session created: $session_name"
    tmux new-session -d -s $session_name "bash $SCRIPTDIR/remove_strains_from_seqid2taxid_v2.sh $nodes_segement $nodes_file $seq_tax_id_map_copy $LOGFILE $output_file"

    # Update the start and finish line numbers for the next chunk
    start=$((start + div))
    fin=$((fin + div))

done
