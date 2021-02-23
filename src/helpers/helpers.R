library(data.table)
library(tidyverse)
library(forcats)
library(coriell)



# process REdiscoverTE output ---------------------------------------------
process_quant_file <- function(quant_file) {
  # set up REdiscover annotations
  gencove_annot <- as.data.table(readRDS("/mnt/data/gdata/human/REdiscoverTE_hg38/REdiscoverTE_hg38_GFP/GENCODE.V26.Basic_Gene_Annotation_md5_GFP.rds"))
  rmsk_annot <- as.data.table(readRDS("/mnt/data/gdata/human/REdiscoverTE_hg38/rmsk_annotation.RDS"))
  rmsk_intergenic <- rmsk_annot[selected_feature == "intergenic"]
  rmsk_intron <- rmsk_annot[selected_feature == "intron"]
  rmsk_exon <- rmsk_exon <- rmsk_annot[selected_feature == "exon"]
  
  DT <- data.table::fread(quant_file, sep = "\t")
  DT <- DT[NumReads > 0]
  
  transcript_counts <- data.table::merge.data.table(
    x = DT,
    y = gencove_annot,
    by.x = "Name",
    by.y = "md5"
  )
  
  repElem_counts <- data.table::merge.data.table(
    x = DT,
    y = rmsk_annot,
    by.x = "Name",
    by.y = "md5"
  )
  
  intergenicRE_counts <- data.table::merge.data.table(
    x = DT,
    y = rmsk_intergenic,
    by.x = "Name",
    by.y = "md5"
  )
  
  intronicRE_counts <- data.table::merge.data.table(
    x = DT,
    y = rmsk_intron,
    by.x = "Name",
    by.y = "md5"
  )
  
  exonicRE_counts <- data.table::merge.data.table(
    x = DT,
    y = rmsk_exon,
    by.x = "Name",
    by.y = "md5"
  )
  
  # summarize counts -----------------------------------------------------------
  gene_counts <- transcript_counts[,
                                   .(count = sum(NumReads, na.rm = TRUE)),
                                   by = .(symbol)
                                   ]
  
  repElem_counts <- repElem_counts[,
                                   .(repElem = paste(repClass, repFamily, repName, sep = "."),
                                     NumReads)][,
                                                .(count = sum(NumReads, na.rm = TRUE)), 
                                                by = .(repElem)]
  
  intergenicRE_counts <- intergenicRE_counts[,
                                             .(repElem = paste(repClass, repFamily, repName, sep = "."),
                                               NumReads)][,
                                                          .(count = sum(NumReads, na.rm = TRUE)), 
                                                          by = .(repElem)]
  
  intronicRE_counts <- intronicRE_counts[,
                                         .(repElem = paste(repClass, repFamily, repName, sep = "."),
                                           NumReads)][,
                                                      .(count = sum(NumReads, na.rm = TRUE)), 
                                                      by = .(repElem)]
  
  exonicRE_counts <- exonicRE_counts[,
                                     .(repElem = paste(repClass, repFamily, repName, sep = "."),
                                       NumReads)][,
                                                  .(count = sum(NumReads, na.rm = TRUE)), 
                                                  by = .(repElem)]
  
  list(
    "gene_counts" = gene_counts,
    "allRE_counts" = repElem_counts,
    "intergenicRE_counts" = intergenicRE_counts,
    "intronicRE_counts" = intronicRE_counts,
    "exonicRE_counts" = exonicRE_counts
  )
}


# perform differential expression on all contrasts ------------------------
get_de_results <- function(contrast_name, 
                           plot_title, 
                           glm_fit, 
                           contrast_matrix, 
                           fdr = 0.1, 
                           fc = 1.5) {
  
  res_df <- glmTreat(glm_fit, contrast = contrast_matrix[, contrast_name], lfc = log2(fc)) %>% 
    coriell::edger_to_df()
  
  vplot <- coriell::plot_volcano(res_df, fdr = fdr, lfc = log2(fc)) +
    ggtitle(plot_title)
  
  md_plot <- coriell::plot_md(res_df, fdr = fdr, lfc = log2(fc)) + 
    ggtitle(plot_title)
  
  list("table" = res_df, "vplot" = vplot, "mdplot" = md_plot)
}


# create dotplots of RE expression ---------------------------------------------
plot_fam_counts = function(df, con, n_dots = 25) {
  plot_df <- df %>% 
    mutate(direction = case_when(FDR < sig & logFC > 0 ~ "up",
                                 FDR < sig & logFC < 0 ~ "down",
                                 TRUE ~ "non-de")) %>% 
    filter(contrast == con) %>% 
    separate(feature_id, into = c("class", "family", "subfamily"), sep = "\\.") %>% 
    group_by(family, direction) %>% 
    summarize(count = n(),
              mean_lfc = mean(logFC),
              .groups = "drop") %>% 
    filter(direction != "non-de") %>%
    mutate(count = if_else(direction == "down", -count, count)) %>% 
    slice_max(order_by = abs(count), n = n_dots)
  
  if (nrow(plot_df) == 0) return(NULL)
  
  break_values <- pretty(plot_df$count)
  
  plot_df %>%
    ggplot(
      aes(
        x = fct_reorder(family, count, function(x) sum(abs(x))),
        y = count,
        color = direction,
        size = abs(mean_lfc)
        )
      ) +
    geom_hline(yintercept = 0, linetype = 2) +
    geom_point() +
    coord_flip() +
    scale_color_manual(values = c("up" = "firebrick", "down" = "steelblue")) +
    scale_y_continuous(breaks = break_values, labels = abs(break_values)) +
    theme_light() +
    labs(title = con,
         x = NULL,
         y = "Count",
         size = "|mean(logFC)|",
         color = "Differential Expression")
}
