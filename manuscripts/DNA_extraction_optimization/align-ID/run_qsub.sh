#!/bin/bash

# module load nextflow
# module load Python/3.11.1

SCRIPTDIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

DATETIME=$(date '+%Y-%m-%d_%H-%M-%S')
EMAIL=$USER"@cdc.gov"

workflow=$1
config=$2

if [ $workflow == "short" ]; then

    qsub \
        -N Align-ID_short_$USER'_'$DATETIME \
        -M $EMAIL \
        -m abe \
        -q all.q \
        -pe smp 8 \
        -l h_vmem=40G \
        -cwd \
        -o $SCRIPTDIR/qsub_logs \
        -e $SCRIPTDIR/qsub_logs/Align-ID_short_$USER'_'$DATETIME'_error.log' \
        qsub_short.sh $config

elif [ $workflow == "long" ]; then

    qsub \
        -N Align-ID_long_$USER'_'$DATETIME \
        -M $EMAIL \
        -m abe \
        -q all.q \
        -pe smp 8 \
        -l h_vmem=40G \
        -cwd \
        -o $SCRIPTDIR/qsub_logs \
        -e $SCRIPTDIR/qsub_logs/Align-ID_long_$USER'_'$DATETIME'_error.log' \
        qsub_long.sh $config
    
else

    echo "ERROR: Invalid [workflow] parameter: $workflow; must be set to either 'long' or 'short'"

fi
