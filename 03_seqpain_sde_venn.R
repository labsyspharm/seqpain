## Sharon Grim
## 2023

## 03_seqpain_sde_venn.R

##This script makes Venn diagrams showing the distribution of 
##differentially expressed genes among treatments and timepoints.
##All scripts assume we are working in a local directory called "seqpain"
##and that data files are within "/seqpain/data/"
##and that figures will be output in "/seqpain/figures/"

library(tidyverse)
library(BiocParallel)
library(RVenn)
library(ggVennDiagram)
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


F_venn_sde <- function(sde_list){
  vars <- sde_list %>% names
  myBreaks <- map(sde_list,
                  ~ (.x %>%
                       purrr::flatten() %>%
                       sapply(., length) %>%
                       unlist %>%
                       max(.) %>%
                       signif(., digits = 2)) 
  ) %>%
    setNames(., vars) %>%
    map(., ~
          c(seq(from = 0, to = .,
                length.out = 101)))
  legendBreaks <- map(myBreaks,
                      ~ c(seq(from = min(.x), to = max(.x), by = ceiling(length(.x)/3)))) %>%
    setNames(., vars)
  list(myBreaks, legendBreaks) %>%
    setNames(., c("myBreaks", "legendBreaks"))
}

F_venn_legend <- function(venn_df, counttype){
  vars <- venn_df %>% names
  myBreaks <- venn_breaks_list[["myBreaks"]] %>% map(., max) %>% unlist()
  
  temp_df <- imap(venn_df,
                  ~ .x %>%
                    as_tibble %>%
                    dplyr::select(venn_name, var2) %>% 
                    pivot_longer(cols = -c(venn_name),
                                 names_to = "metric",
                                 values_to = "variable") %>%
                    dplyr::select(-metric) %>%
                    distinct(.keep_all = TRUE) %>%
                    collect %>%
                    droplevels)
  
  
  map(temp_df %>% names,
      ~list(
        ((ggplot(data = temp_df[[.x]], aes(x = 1, y = 1, color = variable)) +
            geom_blank() +
            geom_point() +
            theme_void() +
            scale_color_manual(values = myColors,
                               breaks = allvars,
                               drop = TRUE,
                               guide = "legend") +
            guides(color = guide_legend(order = 2, title = "Variable", direction = "vertical",
                                        override.aes = list(shape = 21, size = 3,
                                                            stroke = 2)))) %>%
           g_legend()),
        ((ggplot(data = data.frame(gene_count = c(0, myBreaks[[.x]])),
                 aes(x = 1, y = 1, fill = gene_count)) +
            geom_blank() +
            geom_point() +
            theme_void() +
            scale_fill_viridis(discrete = FALSE, option = "mako",
                               limits = c(0, myBreaks[[.x]]),
                               name = {{counttype}}, guide = "colorbar")) %>%
           g_legend())
      ) %>%
        setNames(., c("p1", "p2"))
  ) %>%
    setNames(., vars)
}

patch_venn2 <- function(vlist, colors, myBreaks){
  myColors <- colors
  patch_plots <- map(vlist %>% names,
                     ~vlist[[.x]] %>%
                       RVenn::Venn(.) %>%
                       ggVennDiagram::process_data(.)) %>%
    setNames(., vlist %>% names)
  
  patch_plots <- 
    map(patch_plots %>% names,
        ~(ggplot() +
            geom_sf(aes(fill = count), data = ggVennDiagram::venn_region(patch_plots[[.x]]), show.legend = FALSE) +
            geom_sf(aes(color = name), linewidth = 5, data = ggVennDiagram::venn_setedge(patch_plots[[.x]]), show.legend = FALSE) +
            # geom_sf_text(aes(label = name), size = 3, fontface = "bold", data = ggVennDiagram::venn_setlabel(patch_plots[[.x]])) +
            geom_sf_label(aes(label = count), alpha = 0.5, data = ggVennDiagram::venn_region(patch_plots[[.x]])) +
            scale_fill_viridis(discrete = FALSE, option = "mako",
                               limits = c(0, max(myBreaks)),
                               name = "Gene count") +
            scale_color_manual(values = myColors,
                               breaks = allvars,
                               drop = TRUE,
                               guide = "legend") +
            theme(legend.position = "bottom") +
            labs(title = .x)) +
          theme_void())
  purrr::reduce(patch_plots, `%+%`)
}



# Import DESeq results ----------------------------------------------------

dir <- getwd()

res_nociceptor_log2FC_treatment_list <- read_rds(file = paste0(dir, "/data/res_nociceptor_log2FC_treatment_list.rds"))
res_nociceptor_sig_list_renamed_tpm <- read_rds(file = paste0(dir, "/data/res_nociceptor_sig_list_renamed_tpm.rds"))

res_brca_log2FC_treatment_list <- read_rds(file = paste0(dir, "/data/res_brca_log2FC_treatment_list.rds"))
res_brca_sig_list_renamed_tpm <- read_rds(file = paste0(dir, "/data/res_brca_sig_list_renamed_tpm.rds"))




# Generate Venn datasets --------------------------------------------------

keep <- c("2A", "2K", "2L", "2M")
#nociceptor venn
{
  ##Make the df of venn lists
  temp_df <- res_nociceptor_log2FC_treatment_list[["sde_deseq"]] %>%
    left_join(., res_nociceptor_sig_list_renamed_tpm[["renamer"]] %>%
                dplyr::filter(., grepl(paste(keep, collapse="|"), second, perl = TRUE)) %>%
                droplevels %>%
                dplyr::select(label, cell_line, treatment, time) %>%
                dplyr::mutate(across(.cols = everything(), as.factor))) %>%
    collect %>%
    droplevels
  
  venn_list <- bind_rows(
    (temp_df %>% #looking at comparisons between timepoints  within same treatment
       dplyr::filter(if_any(c(cell_line, treatment), ~ !grepl("_vs_", .x))) %>%
       dplyr::filter(., grepl("_vs_", time)) %>%
       collect %>%
       droplevels %>%
       dplyr::mutate(contrast = "time",
                     constant = treatment,
                     var1 = time)),
    (temp_df %>% #looking at comparisons between treatments within same timepoints
       dplyr::filter(if_any(c(cell_line, time), ~ !grepl("_vs_", .x))) %>%
       dplyr::filter(., grepl("_vs_", treatment)) %>%
       collect %>%
       droplevels %>%
       dplyr::mutate(contrast = "treatment",
                     var1 = treatment,
                     constant = time))) %>%
    collect %>%
    droplevels %>%
    split(., f = .$contrast) %>%
    map(., droplevels) %>%
    map(., ~.x %>%
          split(., f = .$constant) %>%
          map(., ~.x %>%
                split(., f = .$var1) %>%
                map(as_tibble) %>%
                map(., dplyr::select, gene_symbol, var1,
                    contrast, constant) %>%
                map(droplevels) %>%
                map(., ~.x %>%
                      dplyr::mutate(across(everything(), as.character))))
    )
  
  temp_list <- map(venn_list %>% names,
                   ~ setNames(venn_list[[.x]], paste0(names(venn_list[[.x]]), "_", .x))) %>%
    setNames(., names(venn_list))
  

  venn_df <- map(temp_list,
                 ~imap_dfr(., ~ .x %>%
                             bind_rows %>%
                             as_tibble %>%
                             dplyr::mutate(venn_name = .y) %>%
                             dplyr::select(venn_name, contrast, var1, constant, gene_symbol) %>%
                             pivot_longer(, cols = c(gene_symbol),
                                          names_to = "metric",
                                          values_to = "gene_symbol") %>%
                             dplyr::select(-metric) %>%
                             dplyr::mutate(across(everything(), factor)) %>%
                             collect %>%
                             droplevels) %>%
                   split(., f = .$constant))
  
  venn_df <- map(venn_df,
                 ~map(., ~.x %>%
                        droplevels))
  
  ##Make the lists of venn lists
  venn_list_plot <- map(venn_df,
                        ~map(., ~ .x %>%
                               split(., f = .$var1) %>%
                               map(., ~.x %>%
                                     dplyr::select(gene_symbol) %>%
                                     droplevels %>%
                                     dplyr::mutate(across(everything(), as.character)) %>%
                                     unlist) %>%
                               Filter(function(k) length(k) > 0, .) #https://stackoverflow.com/questions/15292904/remove-empty-list-from-a-list-of-lists
                        ) 
  ) %>%
    map(., ~ .x %>%
          Filter(function(k) length(k) > 1, .)) #have to remove contrasts in which var2 has only 1 factor
  
  ##Make the figure legends consistent between venns of the same celltype
  
  myBreaks <- purrr::flatten(venn_list_plot) %>%
    purrr::flatten() %>%
    sapply(., length) %>%
    unlist %>%
    max %>%
    signif(., digits = 2)
  myBreaks <- c(seq(from = 0, to = myBreaks,
                    by = (myBreaks / 100)))
  legendBreaks <- c(seq(from = min(myBreaks), to = max(myBreaks), by = ceiling(length(myBreaks)/3)))
  
  allvars <- intersect((res_nociceptor_sig_list_renamed_tpm[["renamer"]] %>%
                          dplyr::select(treatment, time) %>%
                          bind_rows %>%
                          distinct(.keep_all = TRUE) %>%
                          unlist %>%
                          as.character %>%
                          unique),
                       (imap_dfr(venn_df %>%
                                   purrr::flatten(),
                                 ~ .x %>%
                                   as_tibble %>%
                                   dplyr::select(var1, constant) %>% 
                                   bind_rows %>% 
                                   distinct(.keep_all = TRUE)) %>%
                          unlist %>%
                          as.character %>%
                          unique))
  
  myColors <- viridisLite::turbo(n = length(allvars)) %>%
    setNames(., allvars[rank(allvars)])
  
  temp_df <- imap(venn_df %>% purrr::flatten(),
                  ~ .x %>%
                    as_tibble %>%
                    dplyr::select(contrast, var1) %>% 
                    pivot_longer(cols = -c(contrast),
                                 names_to = "metric",
                                 values_to = "variable") %>%
                    dplyr::select(-metric) %>%
                    distinct(.keep_all = TRUE) %>%
                    collect %>%
                    droplevels) %>%
    bind_rows %>%
    distinct(.keep_all = TRUE) %>%
    split(., f = .$contrast) %>%
    map(., droplevels)
  
  
  venn_legend <- map(temp_df %>% names,
                     ~list(
                       ((ggplot(data = temp_df[[.x]], aes(x = 1, y = 1, color = variable)) +
                           geom_blank() +
                           geom_point() +
                           theme_void() +
                           scale_color_manual(values = myColors,
                                              breaks = allvars,
                                              drop = TRUE,
                                              guide = "legend") +
                           guides(color = guide_legend(order = 2, title = "Variable", direction = "vertical",
                                                       override.aes = list(shape = 21, size = 3,
                                                                           stroke = 2)))) %>%
                          g_legend()),
                       ((ggplot(data = data.frame(gene_count = c(0, max(myBreaks))),
                                aes(x = 1, y = 1, fill = gene_count)) +
                           geom_blank() +
                           geom_point() +
                           theme_void() +
                           scale_fill_viridis(discrete = FALSE, option = "mako",
                                              limits = c(0, max(myBreaks)),
                                              name = "Gene count", guide = "colorbar")) %>%
                          g_legend())
                     ) %>%
                       setNames(., c("p1", "p2"))
  ) %>%
    setNames(., temp_df %>% names)
  

  ##Save the output
  venn_sde_noci <- list(venn_list_plot, venn_df, myBreaks, myColors, venn_legend) %>%
    setNames(., c("venn_list_plot", "venn_df", "myBreaks", "myColors", "venn_legend"))
  
}
##venn_sde_noci

#brca venn
{
  ##Make the df of venn lists
  temp_df <- res_brca_log2FC_treatment_list[["sde_deseq"]] %>%
    left_join(., res_brca_sig_list_renamed_tpm[["renamer"]] %>%
                dplyr::filter(., grepl(paste(keep, collapse="|"), second, perl = TRUE)) %>%
                droplevels %>%
                dplyr::select(label, cell_line, treatment, time) %>%
                dplyr::mutate(across(.cols = everything(), as.factor))) %>%
    collect %>%
    droplevels
  
  venn_list <- bind_rows(
    (temp_df %>% #looking at comparisons between timepoints and treatments within same cell line
       dplyr::filter(., !grepl("_vs_", cell_line)) %>%
       dplyr::filter(., !grepl("_vs_", time)) %>%
       dplyr::filter(., grepl("_vs_", treatment)) %>%
       # dplyr::filter(if_any(c(treatment, time), ~ grepl("_vs_", .x))) %>%
       collect %>%
       droplevels %>%
       dplyr::mutate(contrast = "treatment_time",
                     var1 = treatment,
                     var2 = time,
                     constant = cell_line)),
    (temp_df %>% #looking at comparisons between treatments and cell lines within timepoints
       # dplyr::filter(., grepl("_vs_", cell_line)) %>%
       dplyr::filter(., !grepl("_vs_", time)) %>%
       # dplyr::filter(., grepl("_vs_", treatment)) %>%
       dplyr::filter(if_any(c(cell_line, treatment), ~ grepl("_vs_", .x))) %>%
       collect %>%
       droplevels %>%
       dplyr::mutate(contrast = "cellline_treatment",
                     var1 = cell_line,
                     var2 = treatment,
                     constant = time))) %>%
    collect %>%
    droplevels %>%
    split(., f = .$contrast) %>%
    map(., droplevels) %>%
    map(., ~.x %>%
          split(., f = .$constant) %>%
          map(., ~.x %>%
                split(., f = .$var1) %>%
                map(as_tibble) %>%
                map(., dplyr::select, gene_symbol,
                    var1, var2,
                    contrast, constant) %>%
                map(droplevels) %>%
                map(., ~.x %>%
                      dplyr::mutate(across(everything(), as.character))))
    )
  
  temp_list <- map(venn_list %>% names,
                   ~ setNames(venn_list[[.x]], paste0(names(venn_list[[.x]]), "_", .x)))
  temp_list <- purrr::flatten(temp_list)
  temp_list <- map(temp_list %>% names,
                   ~ setNames(temp_list[[.x]], paste0(names(temp_list[[.x]]), "_", .x))) %>%
    setNames(., names(temp_list))
  
  venn_df <- temp_list %>%
    imap_dfr(., ~ .x %>%
               bind_rows %>%
               as_tibble %>%
               dplyr::mutate(venn_name = .y) %>%
               dplyr::select(venn_name, contrast, var1, var2, constant, gene_symbol) %>%
               pivot_longer(, cols = c(gene_symbol),
                            names_to = "metric",
                            values_to = "gene_symbol") %>%
               dplyr::select(-metric) %>%
               dplyr::mutate(across(everything(), factor)) %>%
               collect %>%
               droplevels) %>%
    split(., f = .$venn_name)
  
  ##Make the lists of venn lists
  
  venn_list_plot <- venn_df %>%
    map(., ~.x %>%
          split(., f = .$var1) %>%
          map(., ~.x %>%
                split(., f = .$var2) %>%
                map(., ~.x %>%
                      dplyr::select(gene_symbol) %>%
                      droplevels %>%
                      dplyr::mutate(across(everything(), as.character)) %>%
                      unlist) %>%
                Filter(function(k) length(k) > 0, .) #https://stackoverflow.com/questions/15292904/remove-empty-list-from-a-list-of-lists
          )
    ) %>%
    map(., ~ .x %>%
          Filter(function(k) length(k) > 1, .)) #have to remove contrasts in which var2 has only 1 factor
  
  ##Make the figure legends consistent between venns of the same celltype
  venn_breaks_list <- F_venn_sde(venn_list_plot)
  
  venn_df <- map(venn_df,
                 droplevels)
  
  allvars <- intersect((res_brca_sig_list_renamed_tpm[["renamer"]] %>%
                          dplyr::select(cell_line, treatment, time) %>%
                          bind_rows %>%
                          distinct(.keep_all = TRUE) %>%
                          unlist %>%
                          as.character %>%
                          unique),
                       (imap_dfr(venn_df,
                                 ~ .x %>%
                                   as_tibble %>%
                                   dplyr::select(var1, var2, constant) %>% 
                                   bind_rows %>% 
                                   distinct(.keep_all = TRUE)) %>%
                          unlist %>%
                          as.character %>%
                          unique))
  
  myColors <- viridisLite::turbo(n = length(allvars)) %>%
    setNames(., allvars[rank(allvars)])

  venn_legend <- F_venn_legend(venn_df, "Gene counts")
  
  
  ##Save the output
  
  venn_sde_brca <- list(venn_list_plot, venn_df, venn_breaks_list, venn_legend, allvars, myColors) %>%
    setNames(., c("venn_list_plot", "venn_df", "venn_breaks_list", "venn_legend", "allvars", "myColors"))
  
}
##venn_sde_brca

##if you want:
##what genes are shared in these venns?
# {
# venn_sde_noci_lists <- map(venn_sde_noci[["venn_list_plot"]] %>% names,
#                            ~map(venn_sde_noci[["venn_list_plot"]][[.x]],
#                                 ~RVenn::Venn(.))) %>%
#   setNames(., venn_sde_noci[["venn_list_plot"]] %>% names)
# 
# 
# venn_sde_brca_lists <- map(venn_sde_brca[["venn_list_plot"]] %>% names,
#                            ~map(venn_sde_brca[["venn_list_plot"]][[.x]],
#                                 ~RVenn::Venn(.))) %>%
#   setNames(., venn_sde_brca[["venn_list_plot"]] %>% names)
# }


write_rds(venn_sde_noci, paste0(dir, "/data/venn_sde_noci.rds"), compress = "gz")
write_rds(venn_sde_brca, paste0(dir, "/data/venn_sde_brca.rds"), compress = "gz")


# Make Venn figures -------------------------------------------------------

#within treatments, and within timepoints:

map(venn_sde_noci[["venn_list_plot"]] %>% names,
                 ~(((patch_venn2(venn_sde_noci[["venn_list_plot"]][[.x]],
                                 venn_sde_noci[["myColors"]],
                                 venn_sde_noci[["myBreaks"]]) +
                       theme(legend.position = "none") +
                       plot_layout(nrow = 2, ncol = 2)) | (wrap_elements(venn_sde_noci[["venn_legend"]][[.x]][["p2"]], ignore_tag = TRUE) / wrap_elements(venn_sde_noci[["venn_legend"]][[.x]][["p1"]], ignore_tag = TRUE))) +
                     plot_layout(ncol = 2, widths = c(3, 1)) +
                     plot_annotation(tag_levels = "A",
                                     title = paste0("Comparing within nociceptor samples")))) %>%
  setNames(., venn_sde_noci[["venn_list_plot"]] %>% names) %>%
  map2(., names(.),
       ~ggsave(paste0(dir, "/figures/", Sys.Date(), "-venn_SDE_noci-", .y, ".png"),
               .x,
               width = 10, height = 7, units = "in"))
###for Brca

map(venn_sde_brca[["venn_list_plot"]] %>% names,
                 ~((patch_venn2(venn_sde_brca[["venn_list_plot"]][[.x]],
                                venn_sde_brca[["myColors"]],
                                venn_sde_brca[["venn_breaks_list"]][["myBreaks"]][[.x]]) +
                      theme(legend.position = "none") +
                      plot_layout(ncol = 2)) | (wrap_elements(venn_sde_brca[["venn_legend"]][[.x]][["p2"]], ignore_tag = TRUE) / wrap_elements(venn_sde_brca[["venn_legend"]][[.x]][["p1"]], ignore_tag = TRUE))) +
                   plot_layout(ncol = 2, widths = c(3, 1)) +
                   plot_annotation(tag_levels = "A",
                                   title = paste0("Comparing within ", paste0(.x)))) %>%
  setNames(., venn_sde_brca[["venn_list_plot"]] %>% names) %>%
  map2(., names(.),
     ~ggsave(paste0(dir, "/figures/", Sys.Date(), "-venn_SDE_brca-", .y, ".png"),
             .x,
             width = 10, height = 7, units = "in"))



