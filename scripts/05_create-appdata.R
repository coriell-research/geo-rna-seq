# Create combined metadata SummarizedExperiment objects
#
# This script will import all SE files from each experiment and extract the DE
# results into a single data.table. Individual columns of the data.table are then
# cast wider into separate matrices, one for each measurement. Metadata information
# annotated in Excel is read in and the metadata and measurement matrices are
# then used to create a single Summarized Experiment of the DE results.
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
de[,  `:=`(id = str_c(BioProject, contrast, sep = "."),
           feature_type = fifelse(feature_id %like% ".*\\..*\\..*", "TE", "Gene"))]

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
  message("Would you like to continue with downstream processing?")
  ans <- menu(c("yes", "no"), title = "Continue Processing?")
  if (ans != 1) stop()
}

if (length(missing_data) > 0) {
  message("There are ", length(missing_data), " IDs annotated without data!")
  message("Writing these IDs to: appdata/missing-data.tsv")
  fwrite(data.table(id = missing_data), here("appdata", "missing-data.tsv"), sep = "\t")
  message("Would you like to continue with downstream processing?")
  ans <- menu(c("yes", "no"), title = "Continue Processing?")
  if (ans != 1) stop()
}

message("Inner joining metadata onto differential expression data...")
metadata <- metadata[drugs, on = "drug", nomatch = 0L]
metadata <- metadata[cells, on = "cell_line", nomatch = 0L]

# Shape into matrices ----------------------------------------------------------
message("Creating assay data from differential expression results...")
lfc <- as.matrix(
  dcast(
    data = de[, .(id, feature_id, logFC)], 
    formula = feature_id ~ id,
    value.var = "logFC",
    fill = 0.0
    ),
  rownames = "feature_id"
)

fdr <- as.matrix(
  dcast(
    data = de[, .(id, feature_id, adj.P.Val)], 
    formula = feature_id ~ id,
    value.var = "adj.P.Val",
    fill = NA
  ),
  rownames = "feature_id"
)

lcpm <- as.matrix(
  dcast(
    data = de[, .(id, feature_id, AveExpr)], 
    formula = feature_id ~ id,
    value.var = "AveExpr",
    fill = NA
  ),
  rownames = "feature_id"
)

tstat <- as.matrix(
  dcast(
    data = de[, .(id, feature_id, t)], 
    formula = feature_id ~ id,
    value.var = "t",
    fill = NA
  ),
  rownames = "feature_id"
)

stderr <- as.matrix(
  dcast(
    data = de[, .(id, feature_id, SE)], 
    formula = feature_id ~ id,
    value.var = "SE",
    fill = NA
  ),
  rownames = "feature_id"
)

# Create SummarizedExperiment object of DE results -----------------------------
setDF(metadata, rownames = metadata$id)
keep <- intersect(colnames(lfc), rownames(metadata))
metadata <- metadata[keep, ]
lfc <- lfc[, keep]
fdr <- fdr[, keep]
tstat <- tstat[, keep]
lcpm <- lcpm[, keep]
stderr <- stderr[, keep]
stopifnot("rownames of metadata do not match colnames of matrices!" = all(colnames(lfc) == rownames(metadata)))

message("Creating final SummarizedExperiment object from differential expression results...")
se <- SummarizedExperiment(
  assays = list(lfc = lfc, 
                fdr = fdr,
                stat = tstat,
                lcpm = lcpm,
                stderr = stderr),
  colData = metadata
)
rowData(se)$feature_type <- fifelse(rownames(se) %like% ".*\\..*\\..*", "TE", "Gene")
saveRDS(se, here("appdata", "se.rds"), compress = FALSE)

# Select input data ---------------------------------------------------------------------------
# create data used by the shiny app for populating select inputs
select_inputs <- lapply(metadata, \(x) sort(unique(x)))
saveRDS(select_inputs, here("appdata", "select-inputs.rds"))
message("Appdata creation complete.")
message("Done.")
