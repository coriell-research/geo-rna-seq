#!/usr/bin/Rscript
# Create combined metadata SummarizedExperiment objects
#
# This script will import all SE files from each experiment and extract the DE
# results into a single SummarizedExperiment object used as the backend database
# for the Shiny app. Metadata information, annotated manually in Excel, is read
# in and aligned with columns of the imported DE data in order to create a 
# finalized SE object.
# ------------------------------------------------------------------------------
message("Loading libraries...")
suppressPackageStartupMessages(library(here))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(SummarizedExperiment))
suppressPackageStartupMessages(library(data.table))
N_CORES <- 12

message("Getting all processed SummarizedExperiment files...")
se_files <- list.files(
  path = here("bioprojects"),
  pattern = "se.rds",
  recursive = TRUE,
  full.names = TRUE
)
names(se_files) <- str_extract(se_files, "PRJNA[0-9]+")
message("Found ", length(se_files), " results.")

# Extract DE results from SE objects and bind into single DT -------------------
extract_de <- function(fpath) {
  se <- readRDS(fpath)
  metadata(se)[["de"]]
}

message("Extracting DE data from SummarizedExperiments...")
de <- rbindlist(
  parallel::mclapply(se_files, extract_de, mc.cores = N_CORES), 
  idcol = "BioProject"
  )

# Add new column as unique identifier "id"
de[, id := str_c(BioProject, contrast, sep = ".")]

# Read in previous contrast-level metadata files -------------------------------
message("Reading in contrast-level metadata...")
metafile <- here("appdata", "contrast-metadata - metadata.tsv")
drugfile <- here("appdata", "contrast-metadata - drugs.tsv")
cellfile <- here("appdata", "contrast-metadata - cells.tsv")

if (any(!file.exists(c(metafile, cellfile, drugfile)))) {
  message("One of: ", metafile, drugfile, " and ", cellfile, "do not exist!")
  stop()
}
metadata <- fread(metafile)
drugs <- fread(drugfile)
cells <- fread(cellfile)

# Which IDs are annotated in the metadata?
analyzed_ids <- de[, unique(id)]
annotated_ids <- metadata[, unique(id)]
missing_annotation <- setdiff(analyzed_ids, annotated_ids)
missing_data <- setdiff(annotated_ids, analyzed_ids)

if (length(missing_annotation) > 0) {
  message("There are ", length(missing_annotation), " IDs without an annotation!")
  message("Writing these IDs to: appdata/missing-annotations.tsv")
  fwrite(data.table(id = missing_annotation), here("appdata", "missing-annotations.tsv"), sep = "\t")
}

if (length(missing_data) > 0) {
  message("There are ", length(missing_data), " IDs annotated without data!")
  message("Writing these IDs to: appdata/missing-data.tsv")
  fwrite(data.table(id = missing_data), here("appdata", "missing-data.tsv"), sep = "\t")
}

message("Inner joining metadata onto differential expression data...")
metadata <- metadata[drugs, on = "drug", nomatch = 0L]
metadata <- metadata[cells, on = "cell_line", nomatch = 0L]

# Shape into matrices ----------------------------------------------------------
message("Creating assay data from differential expression results...")

col2assay <- function(df, rows, cols, vals) {
  m <- dcast(df, get(rows) ~ get(cols), value.var = vals, fill = NA)
  as.matrix(m, rownames = "rows")
}

assay_cols <- c("logFC", "CI.L", "CI.R", "AveExpr", "t", "P.Value", "adj.P.Val", "z")
assays <- vector("list", length(assay_cols))

names(assays) <- assay_cols
for (i in seq_along(assay_cols)) {
  assays[[i]] <- col2assay(de, rows = "feature_id", cols = "id", vals = assay_cols[[i]])
}

# Create SummarizedExperiment object of DE results -----------------------------
setDF(metadata, rownames = metadata$id)
keep <- intersect(colnames(assays[[1]]), rownames(metadata))
metadata <- metadata[keep, ]
assays <- lapply(assays, \(x) x[, keep])
stopifnot("rownames of metadata do not match colnames of matrices!" = all(colnames(assays[[1]]) == rownames(metadata)))

message("Creating final SummarizedExperiment object from differential expression results...")
se <- SummarizedExperiment(
  assays = assays,
  colData = metadata
)
rowData(se)$feature_type <- fifelse(rownames(se) %like% ".*\\..*\\..*", "TE", "Gene")
saveRDS(se, here("appdata", "se.rds"), compress = TRUE)

# Select input data ---------------------------------------------------------------------------
# create data used by the shiny app for populating select inputs
select_inputs <- lapply(metadata, \(x) sort(unique(x)))
saveRDS(select_inputs, here("appdata", "select-inputs.rds"))
message("Appdata creation complete.")
message("Done.")
