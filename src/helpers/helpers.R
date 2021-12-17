library(data.table)
library(edgeR)
library(limma)


# Function to retrieve significance results from limma-trend pipeline
get_treat <- function(contrasts, fit) {
  res_dts <- lapply(contrasts, function(contrast) { 
    as.data.table(topTreat(fit, coef = contrast, number = Inf), keep.rownames = "feature_id") })
  names(res_dts) <- contrasts
  
  # Rename some columns to match edgeR output
  res_dts <- lapply(res_dts, function(dt) setnames(dt, old = c("P.Value", "adj.P.Val"), new = c("PValue", "FDR")))

  return(res_dts)
}

# Function to retrieve significance results from glmTreat test
get_glmTreat <- function(contrasts, fit, contrast_matrix, fc) {
  res_dts <- lapply(contrasts, function(contrast) {
    as.data.table(
      topTags(glmTreat(fit, contrast = contrast_matrix[, contrast], lfc = log2(fc)), n = Inf)$table, 
      keep.rownames = "feature_id")
  })
  names(res_dts) <- contrasts
  
  return(res_dts)
}

