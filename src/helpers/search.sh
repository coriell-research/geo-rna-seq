#!/usr/bin/env bash
###
### Search GEO for the given query. Return a flat file with the 
### BioProject Accession, Number of Samples, and Title of Project
###
### Usage:
###   01_search <query> <out_dir>
###
### Options:
###   <query>   QUERY string used to search GEO
###   <out_dir> Directory to save search results file.  
###   -h        Show this message.
# -----------------------------------------------------------------------------
help() {
        sed -rn 's/^### ?//;T;p' "$0"
    }

if [[ $# == 0 ]] || [[ "$1" == "-h" ]]; then
    help
    exit 1
fi

QUERY=$1
OUT_DIR=$2

esearch \
    -db gds \
    -query "$QUERY" | \
efetch \
    -format docsum | \
xtract \
    -pattern DocumentSummary \
    -element BioProject n_samples title > $OUT_DIR/search-results.tsv
