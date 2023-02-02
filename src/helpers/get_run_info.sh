#!/usr/bin/env bash
### Query the SRA database for each Project Accession returned by the search
### script. Return the run-info.csv file for each which contains information
### about the SRA Accessions, treatments, library prep, etc.
###
### Usage:
###   02_get_run_info.sh <search-results.tsv> <out_dir>
###
### Options:
###   <search-results.tsv>   search-results.tsv file from previos step
###   <out_dir>              Desired output directory
###   -h        Show this message.
# -----------------------------------------------------------------------------
set -Eeuo pipefail

help() {
        sed -rn 's/^### ?//;T;p' "$0"
    }

if [[ $# == 0 ]] || [[ "$1" == "-h" ]]; then
    help
    exit 1
fi

SEARCH_RESULTS=$1
OUT_DIR=$2
MIN_SAMPLES=1
MAX_SAMPLES=100

while read -r BIOPROJECT N_SAMPLES TITLE; do
    if [ $N_SAMPLES -gt $MIN_SAMPLES -a $N_SAMPLES -lt $MAX_SAMPLES ] 
    then
        mkdir -p $OUT_DIR/$BIOPROJECT
        sleep 1
        esearch -db sra -query \"$BIOPROJECT\" < /dev/null | efetch -format runinfo > $OUT_DIR/$BIOPROJECT/run-info.csv
    fi
done < $SEARCH_RESULTS

