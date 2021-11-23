# Create combined metadata SummarizeExperiment objects
#
# This script will import all SE files from each experiment and extract the DE
# results into a single data.table. Individual columns of the data.table are then
# cast wider into separate matrices, one for each measurement. Metadata information
# annotated in Excel is read in and the metadata and measurement matrices are
# then used to create a single Summarized Experiment of the DE results.
# ------------------------------------------------------------------------------
suppressPackageStartupMessages(library(here))
suppressPackageStartupMessages(library(SummarizedExperiment))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(DBI))
N_CORES <- 24


# extract gene annotations
annot <- readRDS(here("doc", "gencode.v26.annotation.cds.rds"))
gene_names <- annot$gene_name

# list the SummarizedExperiment files for each experiment
se_files <- list.files(
  path = here("results"),
  pattern = "se(.*)?.rds",
  recursive = TRUE,
  full.names = TRUE
)

# Remove combined datasets if present
se_files <- se_files[!grepl("combined", se_files)]

# Name each by the BioiProject ID
names(se_files) <- regmatches(se_files, regexpr("PRJNA[0-9]+", se_files))

# Extract DE results from SE objects and bind into single DT -------------------
extract_de_results <- function(fpath) {
  se <- readRDS(fpath)
  de_result <- metadata(se)[["de_results"]]
  return(de_result)
}
de_results <- rbindlist(
  parallel::mclapply(se_files, extract_de_results, mc.cores = N_CORES), 
  idcol = "BioProject", 
  fill = TRUE
  )

# Add new column as unique identifier "id"
de_results[,  `:=`(id = paste(BioProject, contrast, sep = "."),
                   feature_type = fifelse(feature_id %chin% gene_names, "gene", "re"),
                   batch_id = NULL)]

# Write all DE results out to a database ---------------------------------------
con <- dbConnect(RSQLite::SQLite(), here("results", "combined", "data-files", "de.sqlite"))
dbWriteTable(con, "de", de_results, overwrite = TRUE)
dbExecute(con, "CREATE INDEX type_idx ON de (feature_type);")
dbExecute(con, "CREATE INDEX id_idx ON de (id);")
dbExecute(con, "CREATE INDEX project_idx ON de (BioProject);")
dbExecute(con, "CREATE INDEX contrast_idx ON de (contrast);")
dbDisconnect(con)

# Read in previous contrast-level metadata files -------------------------------
metadata <- fread(here("doc", "contrast-metadata", "contrast-metadata - metadata.tsv"))
drugs <- fread(here("doc", "contrast-metadata", "contrast-metadata - drugs.tsv"))
cells <- fread(here("doc", "contrast-metadata", "contrast-metadata - cells.tsv"))

# left join drug info onto metadata
metadata <- merge.data.table(
  x = metadata, 
  y = drugs,
  by.x = "drug",
  by.y = "drug",
  all.x = TRUE,
  all.y = FALSE
  )

# left join cell_info onto metadata
metadata <- merge.data.table(
  x = metadata,
  y = cells,
  by.x = "cell_line",
  by.y = "cell_line",
  all.x = TRUE,
  all.y = FALSE
  )

# Write out the IDs that need annotations
fwrite(unique(de_results[!id %chin% metadata$id, .(id, BioProject, contrast)]), 
       file = here("doc", "contrast-metadata", "contrast-annotation-diff.tsv"),
       sep = "\t"
       )

# Shape into matrices ----------------------------------------------------------
## logFC
lfc_mat <- as.matrix(
  dcast(
    data = de_results[, .(id, feature_id, logFC)], 
    formula = feature_id ~ id,
    value.var = "logFC",
    fill = 0L
    ),
  rownames = "feature_id"
)

## unshrunk logFC
ulfc_mat <- as.matrix(
  dcast(
    data = de_results[, .(id, feature_id, unshrunk.logFC)], 
    formula = feature_id ~ id,
    value.var = "unshrunk.logFC",
    fill = 0L
  ),
  rownames = "feature_id"
)

## logCPM
lcpm_mat <- as.matrix(
  dcast(
    data = de_results[, .(id, feature_id, logCPM)], 
    formula = feature_id ~ id,
    value.var = "logCPM",
    fill = 0L
  ),
  rownames = "feature_id"
)

fdr_mat <- as.matrix(
  dcast(
    data = de_results[, .(id, feature_id, FDR)], 
    formula = feature_id ~ id,
    value.var = "FDR",
    fill = 1L
  ),
  rownames = "feature_id"
)

# Create SummarizedExperiment object of DE results -----------------------------
# coerce metadata to data.frame with rownames
setDF(metadata, rownames = metadata$id)

# subset/reorder metadata based on colnames of matrices (any matrix works)
stopifnot("Metadata information missing for some contrasts. Some data may be dropped from final SE object." = all(colnames(lfc_mat) %in% rownames(metadata)))
metadata <- metadata[colnames(lfc_mat), ]

se <- SummarizedExperiment(
  assays = list(lfc = lfc_mat, 
                fdr = fdr_mat,
                ulfc = ulfc_mat,
                lcpm = lcpm_mat),
  colData = metadata
)
rowData(se)$feature_type <- fifelse(rownames(se) %chin% gene_names, "gene", "re")

## REMOVE OUTLIERS BEFORE WRITING OUT ------------------------------------------
outlier_contrasts <- c("PDX_CBAB22694_CD45pos.1percent_VTP50469_chow.28d_vs_control",
                       "Bortezomib_vs_DMSO_24hr.10A",
                       "Bortezomib_vs_DMSO_24hr.10B",
                       "Ciclopirox_vs_DMSO_24hr.12A",
                       "Manumycin_A_vs_DMSO_6hr.12A",
                       "SCH_79797_vs_DMSO_24hr.12A",
                       "GSK461364_vs_DMSO_24hr.10B",
                       "Sorafenib_vs_DMSO_24hr.10B",
                       "TW_37_vs_DMSO_24hr.10B",
                       "Vorinostat_vs_DMSO_24hr.10B"
                       )
se <- se[, !se$contrast %in% outlier_contrasts]

# save object out to be read into interactive visualization
saveRDS(se, here("results", "combined", "rds-files", "se.rds"), compress = FALSE)
