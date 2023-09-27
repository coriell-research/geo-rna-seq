#!/usr/local/bin/Rscript
#
# Create a MultiAssayExperiment object a user supplied annotation file
# and the quant.sf files for each sample processed by the REdiscoverTE pipeline.
#
# NOTE: There are hard-coded paths to annotation files
# ------------------------------------------------------------------------------
suppressMessages(library(optparse, warn.conflicts = FALSE, quietly = TRUE))
suppressMessages(library(data.table, warn.conflicts = FALSE, quietly = TRUE))
suppressMessages(library(MultiAssayExperiment, warn.conflicts = FALSE, quietly = TRUE))
suppressMessages(library(tximport, warn.conflicts = FALSE, quietly = TRUE))

# get commandline arguments
option_list <- list(
  make_option(c("-d", "--quants_dir"),
    type = "character",
    default = NULL,
    help = "Path to quants/ directory containing sub-directories for each SAMPLE",
    metavar = "quants_dir"
  ),
  make_option(c("-f", "--annotation"),
    type = "character",
    default = NULL,
    help = "Path to annotation file. Must at least contain a column called 'Run' and 'group'. Other metadata is optional",
    metavar = "annotation"
  ),
  make_option(c("-o", "--out_dir"),
    type = "character",
    default = ".",
    help = "Location to save the final MultiAssayExperiment object RDS file",
    metavar = "out_dir"
  )
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

message("Reading in annotation files for processing steps...")
tx2gene <- readRDS("/mnt/data/gdata/human/REdiscoverTE_hg38/tx2gene_REdiscoverTE.rds")

message("Creating count matrices from quant files...")
quant_files <- list.files(
  path = opt$quants_dir,
  pattern = "quant.sf",
  recursive = TRUE,
  full.names = TRUE
)
names(quant_files) <- regmatches(quant_files, regexpr("SRR[0-9]+", quant_files))
txi <- tximport(
  files = quant_files,
  type = "salmon",
  countsFromAbundance = "lengthScaledTPM",
  tx2gene = tx2gene,
  importer = function(x) data.table::fread(x, showProgress = FALSE)
)

# Extract individual matrices
intron_mat <- txi$counts[grepl("__intron$", rownames(txi$counts)), ]
intergenic_mat <- txi$counts[grepl("__intergenic", rownames(txi$counts)), ]
exon_mat <- txi$counts[grepl("__exon", rownames(txi$counts)), ]
all_re <- Reduce(union, list(rownames(intron_mat), rownames(intergenic_mat), rownames(exon_mat)))
gene_mat <- txi$counts[!rownames(txi$counts) %chin% all_re, ]

# Combine the intronic and intergenic counts into a single RE matrix
intron_dt <- as.data.table(intron_mat, keep.rownames = "feature_id")
intergenic_dt <- as.data.table(intergenic_mat, keep.rownames = "feature_id")
re_dt <- rbind(intron_dt, intergenic_dt)
re_dt[, feature_id := gsub("__intron$|__intergenic$", "", feature_id)]
re_dt.m <- melt(re_dt, id.vars = "feature_id", variable.name = "Run", value.name = "count")
by_re <- re_dt.m[, .(count = sum(count)), by = .(Run, feature_id)]
re_mat <- as.matrix(dcast(by_re, feature_id ~ Run, value.var = "count", fill = 0.0), rownames = "feature_id")

# Read in the sample metadata
message("Reading in sample annotation information...")
meta_df <- fread(opt$annotation)
stopifnot("Run" %in% colnames(meta_df))
stopifnot("group" %in% colnames(meta_df))

# convert to data.frame for colData
setDF(meta_df, rownames = meta_df$Run)
meta_df <- subset(meta_df, select = -Run)
stopifnot(all(names(quant_files) %in% rownames(meta_df)))

# reorder the metadata to match colnames of matrices
meta_df <- meta_df[names(quant_files), ]

# create MultiAssayExperiment object from each assay and metadata
exp_list <- list(
  "gene counts" = gene_mat,
  "RE counts" = re_mat,
  "exon RE counts" = exon_mat,
  "intron RE counts" = intron_mat,
  "intergenic RE counts" = intergenic_mat
)

message("Creating MultiAssayExperiment object...")
mae <- MultiAssayExperiment(experiments = exp_list, colData = meta_df)

# write mae .rds file to out_dir
message(paste("Writing MultiAssayExperiment object to", file.path(opt$out_dir, "MultiAssayExperiment.rds")))
saveRDS(mae, file = file.path(opt$out_dir, "MultiAssayExperiment.rds"))
message("Done.")
