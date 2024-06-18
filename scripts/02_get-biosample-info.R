#!/usr/bin/Rscript
#
# Query NCBI BioSample for Metadata about each of the samples
#
# ------------------------------------------------------------------------------
suppressPackageStartupMessages(library(optparse, warn.conflicts = FALSE, quietly = TRUE))
suppressPackageStartupMessages(library(data.table, warn.conflicts = FALSE, quietly = TRUE))
suppressPackageStartupMessages(library(rentrez, warn.conflicts = FALSE, quietly = TRUE))
suppressPackageStartupMessages(library(xml2, warn.conflicts = FALSE, quietly = TRUE))

# get commandline arguments
option_list <- list(
  make_option(c("-f", "--run_info"),
              type = "character",
              default = NULL,
              help = "Path to run_info file for this BioProject",
              metavar = "run_info"
  ),
  make_option(c("-o", "--out_dir"),
              type = "character",
              default = ".",
              help = "Location to save the biosample-info and metadata files",
              metavar = "out_dir"
  )
)
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)
outfile <- file.path(opt$out_dir, "biosample-info.tsv")
outfile2 <- file.path(opt$out_dir, "metadata.txt")

# Get unique biosamples from run_info
run_info <- fread(opt$run_info)
biosamples <- unique(run_info[TaxID == "9606" & LibraryStrategy == "RNA-Seq", BioSample])

# Perform a query on ncbi biosample
get_biosample_dt <- function(acc) {
  message(paste("Searching for", acc))
  Sys.sleep(1)
  search_res <- entrez_search(db = "biosample", term = acc)
  search_id <- search_res$ids
  message("Done.")
  message("Getting summary for search ID", search_id)
  summary_list <- entrez_summary("biosample", search_res$ids)
  message("Done.")
  message("Converting results to data.table")
  
  return(as.data.table(summary_list))
}

message("Querying NCBI BioSample for BioSample metadata...")
biosample_dt <- rbindlist(lapply(biosamples, get_biosample_dt), fill = TRUE, use.names = TRUE)

message("Cleaning BioSample Metadata...")
biosample_dt[, c("bio_sample", "sra", "geo") := tstrsplit(identifiers, ";")][,
             `:=`(date = NULL, publicationdate = NULL, modificationdate = NULL,
                 organism = NULL, sourcesample = NULL, identifiers = NULL,
                 infraspecies = NULL, package = NULL, sortkey = NULL, bio_sample = NULL)][,
                `:=`(sra = trimws(gsub("SRA: ", "", sra)), geo = trimws(gsub("GEO: ", "", geo)))]

message("Parsing the BioSample xml attribute strings...")
parse_attrs <- function(col) {
  out <- tryCatch(
    {
      parsed <- xml_find_all(read_xml(col, options = c("RECOVER", "NOERROR")), ".//Attribute")
      val_dt <- data.table(vals = unlist(as_list(parsed)))
      key_dt <- data.table(keys = unlist(lapply(xml_attrs(parsed), function(x) x[[1]])))
      dt <- cbind(key_dt, val_dt)
      return(dt)
    },
    error = function(cond) {
      message("Bad XML")
      message(cond)
      return(data.table(keys = NA, vals = NA))
    }
  )
  return(out)
}
xml_strings <- biosample_dt$sampledata
names(xml_strings) <- biosample_dt$accession
attr_dt <- rbindlist(lapply(xml_strings, parse_attrs), idcol = "BioSample")
attr_dt.w <- dcast(attr_dt, BioSample ~ keys, value.var = "vals")

message("Joining Attributes onto BioSample Metadata...")
biosample_metadata <- merge(
  x = biosample_dt,
  y = attr_dt.w,
  by.x = "accession",
  by.y = "BioSample"
)

# drop xml string
biosample_metadata[, sampledata := NULL]

message(paste("Writing BioSample data to", outfile))
fwrite(x = biosample_metadata, file = outfile, sep = "\t")

# Join biosample data onto run-info
message("Creating metadata from run-info and biosample data...")
metadata <- merge(
  x = run_info,
  y = biosample_metadata,
  by.x = "BioSample",
  by.y = "accession",
  all.x = TRUE,
  all.y = FALSE)
metadata <- metadata[TaxID == "9606" & LibraryStrategy == "RNA-Seq", ][order(BioSample)]

message(paste("Writing metadata to", outfile2))
fwrite(x = metadata, file = outfile2, sep = "\t")

