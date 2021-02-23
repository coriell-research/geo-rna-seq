# GEO Cancer Dataset RNA-Seq Analysis

This directory contains the analysis of RNA-Seq data downloaded from GEO (or 
possibly other sources) processed using the REdiscoverTE pipeline with the goal 
of comparing gene expression across different cell-lines/chemotherapy drugs. 

## data/ 

Contains raw data downloaded from GEO using the `automate-geo-search` scripts. 
Each directory is named by its BIOPROJECT ID and contains gene counts (quants)
for each samples in the BIOPROJECT as well as run-info and additional metadata
downloaded from SRA

See the `~/data/tools/automate-geo-search/` directory for details about scripts
used to download all raw data.

## doc/ 

## results/

Contains results from processing steps. Separate directories for each BIOPROJECT

## src/

Contains source code for each analysis.
