#!/bin/bash
SQLITE3_CONTAINER="https://depot.galaxyproject.org/singularity/sqlite%3A3"

local_nodes_db=$1
nodes=$2
my_tax_ids=$3

bind_dir="/"$(echo $(pwd) | cut -d '/' -f2)


while IFS= read -r id; do #iterating through tax ids listed in the tax id input file

    data=$(singularity exec --bind $bind_dir $SQLITE3_CONTAINER sqlite3 $local_nodes_db "SELECT parent_id FROM TAX_IDS WHERE parent_id = "$id"")


    if [[ -n $data ]]; then

        continue

    fi

    singularity exec --bind $bind_dir $SQLITE3_CONTAINER sqlite3 $local_nodes_db "INSERT INTO TAX_IDS (parent_id, child_ids) VALUES ("$id", '');"

    parent_id="$id"

    tax_array=("$parent_id") #creating array for tax ids
    loop_again=true
    count_other=0

    while $loop_again; do

        count=0

        for i in "${tax_array[@]}"; do

            if [[ $count -eq 0 ]]; then

                tax_array=() #clearing out array to rebuild it with child tax ids for the next for loop iteration
                ((count++))

            fi

            child_ids="$(awk -v id="$i" '$3 == id { print $1 }' $nodes)" #finding child nodes that have the "parent tax id" column set to current tax id

            if [[ ! -z $child_ids ]]; then #checking to see if we found any child nodes

                for x in $child_ids; do

                    if [[ $count_other -eq 0 ]]; then

                        # this should only execute once per specified parent tax id
                        singularity exec --bind $bind_dir $SQLITE3_CONTAINER sqlite3 $local_nodes_db "UPDATE TAX_IDS SET child_ids = "$x" WHERE parent_id = "$parent_id";"
                        ((count_other++))

                        tax_array+=("$x") #appending child tax ids to array to use in
                    
                    else

                        singularity exec --bind $bind_dir $SQLITE3_CONTAINER sqlite3 $local_nodes_db "UPDATE TAX_IDS SET child_ids = child_ids || ',' || "$x" WHERE parent_id = "$parent_id";"

                        tax_array+=("$x") #appending child tax ids to array to use in

                    fi

                done                        

            else 

                continue

            fi

        done

        if [ ${#tax_array[@]} -eq 0 ]; then

            loop_again=false

        fi

    done

done < "$my_tax_ids"

touch complete.txt
