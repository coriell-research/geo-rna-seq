# Combine all results SummarizedExperiment counts and logcounts into a single database
#
#
# -------------------------------------------------------------------------------------------------
suppressPackageStartupMessages(library(here))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(SummarizedExperiment))


# List files for analyzed BioProjects ---------------------------------------------------------

rmd_files <- list.files(
  path = here("src"),
  pattern = "*.Rmd",
  full.names = TRUE
  )
  
# Extract just the project names from the analysis files
projects <- gsub("\\.Rmd", "", basename(rmd_files))

# Remove PRJNA379584 which has multiple batches -- read these in separately
projects <- projects[projects != "PRJNA379584"]

# Create a list of rds-files to import based on analysis scripts
se_files <- here("results", projects, "rds-files", "se.rds")
names(se_files) <- projects

# Get the files for the experiment with multiple batches
se_files2 <- list.files(
  path = here("results", "PRJNA379584", "rds-files"),
  pattern = "se-.*.rds",
  full.names = TRUE
)
names(se_files2) <- gsub("\\.rds", "", gsub("se-", "", basename(se_files2)))
names(se_files2) <- paste0("PRJNA379584_", names(se_files2))

# Combine all into a single vector
se_files_all <- c(se_files, se_files2)

# Read all into a single object
ses <- lapply(se_files_all, readRDS)

