library(data.table)
library(forcats)
library(coriell)
library(tidyverse)


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

# version without the plots
get_de_results2 <- function(contrast_name, 
                            glm_fit, 
                            contrast_matrix, 
                            fc = 1.5) {
  res_df <- glmTreat(
    glm_fit, 
    contrast = contrast_matrix[, contrast_name], 
    lfc = log2(fc)
    ) %>%
    coriell::edger_to_df()
  
  list("table" = res_df)
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
