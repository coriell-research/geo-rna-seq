# GEO Search

The general workflow goes as follows:

1. Find an appropriate study on NCBI and get the BioProject ID (PRJNA[0-9]+)
2. Run `get_run_info_by_PRJNA.sh` to get the run-info.csv file. e.g. `get_run_info_by_PRJNA.sh PRJNA12345 .`
3. Run `get_biosample_info.R` using the run-info.csv as input to get the biosample metadata. e.g. `Rscript get_bio_sample_info.R -f run-info.csv -o .`
  - Add a new column to the metadata.txt output file for group membership of each sample and save as 'annotation.tsv'. The annotation file must contain a 'group' column.
4. Run `download_and_process.py` to download and map all runs in the run-info.csv file. e.g. `python download_and_process.py run-info.csv --threads 8`
5. Run `create_MAE.R` to import the mapped counts into a count matrix for downstream analysis. e.g. `Rscript create_MAE.R -d quants -f annotation.tsv -o . -c 8`
6. Perform DE testing in Rmd file.
  - Copy the template from and older file and change the design matrix. 
7. Import all DE results into combined SummarizedExperiment object
8. Update meta-analysis annotation file on Google Sheets.

