#!/usr/bin/env bash
###
### Download the run-info.csv and metadata.txt file for the given PRJNA ID
###
### Usage:
###   get_run_info_by_PRJNA.sh <PRJ_ID> <out_dir>
###
### Options:
###   <PRJ_ID>   BioProject ID. e.g. PRJNA12345
###   <out_dir>  Directory to save run-info.csv and metadata.txt file.
###   -h         Show this message.
# -----------------------------------------------------------------------------
set -Eeuo pipefail

help() {
        sed -rn 's/^### ?//;T;p' "$0"
    }

if [[ $# == 0 ]] || [[ "$1" == "-h" ]]; then
    help
    exit 1
fi

BIOPROJECT=$1
OUT_DIR=$2

# get the SRA run-info for a single PROJECT ID
esearch -db sra -query \"$BIOPROJECT\" | efetch -format runinfo > $OUT_DIR/run-info.csv

