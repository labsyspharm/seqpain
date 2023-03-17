## Sharon Grim
## 2023

## 05_seqpain_fgsea_hm.R

##This script makes heatmaps of the enrichment scores of pathways
##identified through GSEA.

##All scripts assume we are working in a local directory called "seqpain"
##and that data files are within "/seqpain/data/"
##and that figures will be output in "/seqpain/figures/"

library(tidyverse)
library(BiocParallel)
library(circlize)
library(dendextend)
library(ggplotify)
library(pheatmap)
library(RColorBrewer)
library(viridis)
library(viridisLite)
library(patchwork)


# Custom functions --------------------------------------------------------

g_legend <- function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}


#in this heatmap function, you can choose to either cluster the pathways by:
#gene membership (default)
#enrichment scores (currently commented out)
P_fgsea_hm <- function(fgsealist, renamelist, cutoff){
  metric <- unique(fgsealist$subset)
  namevar <- deparse(substitute(fgsealist)) %>%
    gsub("(^[[:alpha:]]+\\_[[:alpha:]]+).*", "\\1", ., perl = TRUE)
  
  df <- fgsealist %>%
    tidyr::drop_na() %>%
    left_join(., MSigDB_C2_CP_pathways_renamer) %>%
    pivot_wider(., names_from = "renamed",
                id_cols = "systematicName",
                values_from = "NES",
                values_fill = 0) %>% 
    column_to_rownames(var = "systematicName")
  myBreaks <- df %>%
    max(abs(min(.)), abs(max(.))) %>%
    ceiling()
  myBreaks <- c(seq(from = -myBreaks, to = myBreaks,
                    length = 100)) #make the breaks center around 0
  legendBreaks <- c(seq(from = min(myBreaks), to = max(myBreaks), length.out = (max(myBreaks) + 1)))
  
  sample_idx <- renamelist %>%
    dplyr::filter(renamed %in% colnames(df)) %>%
    droplevels %>%
    arrange(., by = cell_type, cell_line, treatment, time) %>%
    dplyr::select(renamed) %>%
    distinct %>%
    simplify %>%
    as.character
  
  annotation_sample <- (renamelist %>%
                          dplyr::select(renamed, cell_line, treatment, time) %>%
                          dplyr::filter(renamed %in% colnames(df)) %>%
                          collect %>%
                          droplevels) %>%
    dplyr::mutate(time = factor(time, levels = c("0h", "1h_vs_0h",
                                                 "1h", "1h_vs_6h", "1h_vs_24h",
                                                 "6h_vs_0h", "6h_vs_1h", "6h", "6h_vs_24h",
                                                 "24h_vs_0h", "24h_vs_1h", "24h_vs_6h", "24h")),
                  treatment = factor(treatment, levels = c("DMSO", "Taxol", "Kw", "KwTaxol",
                                                           "Taxol_vs_DMSO", "Kw_vs_DMSO", "KwTaxol_vs_DMSO",
                                                           "KwTaxol_vs_Taxol", "KwTaxol_vs_Kw", "Taxol_vs_Kw")),
                  cell_line = factor(cell_line, levels = c("iPSCs", "Hs578T", "MDAMB231", "MDAMB231_vs_Hs578T", "MDAMB231_vs_MCF7", 
                                                           "MCF7", "MCF7_vs_Hs578T"))) %>%
    collect %>%
    droplevels %>%
    arrange(., by = cell_line, treatment, time) %>%
    # dplyr::mutate(renamed = fct_relevel(renamed, .$renamed[order(.$renamed, .$cell_line, .$treatment, .$time)])) %>%
    dplyr::mutate(renamed = factor(renamed, levels = unique(renamed))) %>%
    dplyr::rename(sample = renamed) %>%
    collect %>%
    column_to_rownames(var = "sample")
  
  sample_colors <- list(
    cell_line = grey((seq(1:length(levels(annotation_sample$cell_line))) - 0.25) / length(levels(annotation_sample$cell_line)))  %>%
      setNames(., levels(annotation_sample$cell_line)),
    treatment = viridisLite::rocket(n = length(levels(annotation_sample$treatment))) %>%
      setNames(., levels(annotation_sample$treatment)),
    time = viridisLite::viridis(n = length(levels(annotation_sample$time))) %>%
      setNames(., levels(annotation_sample$time))
  )
  
  # dend_pathways <- dist(df, method = "manhattan") %>%
  #   hclust(method = "ward.D2") %>%
  #   as.dendrogram() %>%
  #   ladderize(right = TRUE)
  
  df <- fgsealist %>%
    tidyr::drop_na() %>%
    dplyr::mutate(NES = replace(NES, which(padj > {{cutoff}}), NA)) %>%
    left_join(., MSigDB_C2_CP_pathways_renamer) %>%
    pivot_wider(., names_from = "renamed",
                id_cols = "systematicName",
                values_from = "NES") %>% 
    column_to_rownames(var = "systematicName") %>%
    dplyr::filter(if_any(everything(), ~ !is.na(.)))
  
  #arrange pathways by previous clustering based on gene presence
  # df <- df[match(labels(dend_C2_CP_pathways)[(labels(dend_C2_CP_pathways) %in% rownames(df))], rownames(df)),]
  
  
  #arrange samples by cell_line, treatment, timepoint
  df <- df[match(labels(dend_C2_CP_pathways)[(labels(dend_C2_CP_pathways) %in% rownames(df))], rownames(df)), match(sample_idx[(sample_idx %in% colnames(df))], colnames(df))]
  
  g <- pheatmap(df,
                color = colorRampPalette(pals::ocean.thermal(n = 3))(100),
                border_color = "black",
                scale = "none",
                na_col = "white",
                cluster_rows = FALSE, #if you want to use the clustering of pathways based on gene membership
                # cluster_rows = as.hclust(dend_pathways), #if you want to cluster pathways based on enrichment scores
                annotation_col = annotation_sample, 
                annotation_colors = sample_colors,
                cluster_cols = FALSE, show_colnames = TRUE,
                angle_col = 45, 
                fontsize = 8,
                breaks = myBreaks,
                legend_breaks = legendBreaks,
                main = "", 
                legend_labels = c(as.character(head(legendBreaks, -1)), paste0({metric}, "\n")),
                legend = TRUE,
                silent = TRUE,
                show_rownames = FALSE
  ) 
  legend_list <- grep("legend", g[["gtable"]][["layout"]][["name"]]) 
  legend_list <- map(legend_list, 
                     ~.x %>% 
                       g$gtable$grobs[[.]] %>%
                       as.ggplot()) %>%
    setNames(., g[["gtable"]][["layout"]][["name"]][legend_list])
  legend_list[["annotation_legend"]] <- legend_list[["annotation_legend"]] + 
    theme(plot.margin = unit(c(30,30,5,5), "pt")) + 
    labs(title = paste0("Conditions")) + 
    theme(plot.title = element_text(size = 8, face = "bold"))
  legend_list[["legend"]] <- legend_list[["legend"]] + 
    theme(plot.margin = unit(c(5,30,30,5), "pt")) + 
    labs(title = paste0("NES calculated from ", "\n", {{metric}})) + 
    theme(plot.title = element_text(size = 8, face = "bold", hjust = 0))
  
  g <- pheatmap(df,
                color = colorRampPalette(pals::ocean.thermal(n = 3))(100),
                border_color = "black",
                scale = "none",
                na_col = "white",
                cluster_rows = FALSE,
                annotation_col = annotation_sample, 
                annotation_colors = sample_colors,
                cluster_cols = FALSE, 
                show_colnames = FALSE,
                angle_col = 45, 
                fontsize = 8,
                breaks = myBreaks,
                legend_breaks = legendBreaks,
                main = "", 
                legend_labels = c(as.character(head(legendBreaks, -1)), paste0({metric}, "\n")),
                legend = FALSE,
                annotation_legend = FALSE,
                annotation_names_col = FALSE,
                silent = TRUE,
                show_rownames = FALSE
  ) %>%
    as.ggplot(.) +
    purrr::reduce(legend_list, `/`) +
    plot_layout(widths = c(3,1))
  
}

# Import GSEA results  ----------------------------------------------------

dir <- getwd()

res_nociceptor_fgsea_sig_wide <- read_rds(file = paste0(dir, "/data/res_nociceptor_fgsea_sig_wide.rds"))
res_nociceptor_sig_list_renamed_tpm <- read_rds(file = paste0(dir, "/data/res_nociceptor_sig_list_renamed_tpm.rds"))

res_brca_fgsea_sig_wide <- read_rds(file = paste0(dir, "/data/res_brca_fgsea_sig_wide.rds"))
res_brca_sig_list_renamed_tpm <- read_rds(file = paste0(dir, "/data/res_brca_sig_list_renamed_tpm.rds"))

# Import reference files --------------------------------------------------

MSigDB_C2_CP_pathways_flat <- read_rds(file = paste0(dir, "/data/MSigDB_C2_CP_pathways_flat.rds"))
load(file = paste0(dir, "/data/MSigDB_C2_CP_extra.RData"))
#loads in the following data:
#MSigDB_C2_CP_pathways_v2022_1
#MSigDB_C2_CP_pathways_v2022_genes
#MSigDB_C2_CP_pathways_renamer
#MSigDB_C2_CP_pathways_wide

dend_C2_CP_pathways <- MSigDB_C2_CP_pathways_wide %>%
  dist(., method = "manhattan") %>%
  hclust(method = "ward.D2") %>%
  as.dendrogram() %>%
  ladderize(right = TRUE)
labels(dend_C2_CP_pathways) %>% length

write_rds(dend_C2_CP_pathways, paste0(dir, "/data/dend_C2_CP_pathways.rds"), compress = "gz")


# Generate heat maps: brca ------------------------------------------------------

g1 <- P_fgsea_hm(fgsealist = res_brca_fgsea_sig_wide[["all_log2FC_tpm"]], 
                 renamelist = res_brca_sig_list_renamed_tpm[["renamer"]],
                 cutoff = 1)
g2 <- P_fgsea_hm(fgsealist = res_brca_fgsea_sig_wide[["all_log2FC_tpm"]],
                 renamelist = res_brca_sig_list_renamed_tpm[["renamer"]],
                 cutoff = 0.05)
g3 <- P_fgsea_hm(fgsealist = res_brca_fgsea_sig_wide[["all_log2FC_deseq"]], 
                 renamelist = res_brca_sig_list_renamed_tpm[["renamer"]],
                 cutoff = 1)
g4 <- P_fgsea_hm(fgsealist = res_brca_fgsea_sig_wide[["all_log2FC_deseq"]], 
                 renamelist = res_brca_sig_list_renamed_tpm[["renamer"]],
                 cutoff = 0.05)

g1_legend <- g1[[-1]]
g1[[-1]] <- NULL
g2[[-1]] <- NULL
g3_legend <- g3[[-1]]
g3[[-1]] <- NULL
g4[[-1]] <- NULL

gpatch <- ((g1 | inset_element(g1_legend, align_to = "plot",
                               left = 0.8, bottom = 0, right = 1.3, top = 1, 
                               ignore_tag = TRUE)) | g2 ) / 
  ((g3 | inset_element(g3_legend, align_to = "plot",
                       left = 0.8, bottom = 0, right = 1.3, top = 1, 
                       ignore_tag = TRUE)) | g4) 

gpatch <- gpatch + plot_annotation(
  title = "GSEA",
  caption = paste0("Left: normalized enrichment scores of all pathways", "\n",
                   "Right: only significant pathways"))


ggsave(paste0(dir, "/figures/", Sys.Date(), "-fgsea_hm_brca.png"),
       gpatch,
       width = 16, height = 14, units = "in")


# Generate heatmaps: nociceptor -------------------------------------------

g1 <- P_fgsea_hm(renamelist = res_nociceptor_sig_list_renamed_tpm[["renamer"]],
                 fgsealist = res_nociceptor_fgsea_sig_wide[["all_log2FC_tpm"]],
                 cutoff = 1)
g2 <- P_fgsea_hm(renamelist = res_nociceptor_sig_list_renamed_tpm[["renamer"]],
                 fgsealist = res_nociceptor_fgsea_sig_wide[["all_log2FC_tpm"]],
                 cutoff = 0.05)
g3 <- P_fgsea_hm(renamelist = res_nociceptor_sig_list_renamed_tpm[["renamer"]],
                 fgsealist = res_nociceptor_fgsea_sig_wide[["all_log2FC_deseq"]],
                 cutoff = 1) 
g4 <- P_fgsea_hm(renamelist = res_nociceptor_sig_list_renamed_tpm[["renamer"]],
                 fgsealist = res_nociceptor_fgsea_sig_wide[["all_log2FC_deseq"]],
                 cutoff = 0.05)
g1_legend <- g1[[-1]]
g1[[-1]] <- NULL
g2[[-1]] <- NULL
g3_legend <- g3[[-1]]
g3[[-1]] <- NULL
g4[[-1]] <- NULL

gpatch <- ((g1 | inset_element(g1_legend, align_to = "plot",
                               left = 0.8, bottom = 0.1, right = 1.3, top = 1, 
                               ignore_tag = TRUE)) | g2 ) / 
  ((g3 | inset_element(g3_legend, align_to = "plot",
                       left = 0.8, bottom = 0.1, right = 1.3, top = 1, 
                       ignore_tag = TRUE)) | g4) 
gpatch <- gpatch + plot_annotation(
  title = "GSEA",
  caption = paste0("Left: normalized enrichment scores of all pathways", "\n",
                   "Right: only NES scores of significant pathways are shown"))


ggsave(paste0(dir, "/figures/", Sys.Date(), "-fgsea_hm_noci.png"),
       gpatch,
       width = 16, height = 12, units = "in")


