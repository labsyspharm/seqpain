## Sharon Grim
## 2023

## 04_seqpain_fgsea_venn.R

##This script makes Venn diagrams showing the distribution of 
##significantly up or down-regulated pathways among treatments and timepoints.
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

res_nociceptor_fgsea_sig <- read_rds(file = paste0(dir, "/data/res_nociceptor_fgsea_sig.rds"))
res_nociceptor_sig_list_renamed_tpm <- read_rds(file = paste0(dir, "/data/res_nociceptor_sig_list_renamed_tpm.rds"))

res_brca_fgsea_sig <- read_rds(file = paste0(dir, "/data/res_brca_fgsea_sig.rds"))
res_brca_sig_list_renamed_tpm <- read_rds(file = paste0(dir, "/data/res_brca_sig_list_renamed_tpm.rds"))



# Generate Venn datasets --------------------------------------------------

#nociceptor
{
  temp_df <- map(res_nociceptor_fgsea_sig,
                 ~left_join(.x, res_nociceptor_sig_list_renamed_tpm[["renamer"]] %>%
                              dplyr::select(label, cell_line, treatment, time) %>%
                              collect %>%
                              droplevels %>%
                              dplyr::mutate(across(.cols = everything(), as.factor))) %>% 
                   dplyr::select(renamed, label, cell_line, treatment, time, subset, pathway, padj, NES) %>% 
                   collect %>%
                   droplevels)
  
  #let's keep only "all_log2FC_tpm"
  venn_list <- temp_df[["all_log2FC_tpm"]] %>%
    bind_rows((dplyr::filter(., grepl("_vs_", treatment)) %>%
                 collect %>%
                 droplevels %>%
                 dplyr::mutate(contrast = "longitudinal",
                               var1 = time,
                               var2 = treatment,
                               constant = cell_line)),
              (dplyr::filter(., grepl("_vs_", treatment)) %>% #do you want a treatment-specific patch that has all timepoints for that treatment?
                 collect %>%
                 droplevels %>%
                 dplyr::mutate(contrast = "latitudinal",
                               var1 = treatment,
                               var2 = time,
                               constant = cell_line))) %>%  #timepoint-specific patch that has all treatments within that timepoint
    collect %>%
    droplevels %>%
    split(., f = .$contrast) %>%
    map(., droplevels) %>%
    map(., ~.x %>%
          split(., f = .$var1) %>%
          map(as_tibble) %>%
          map(., dplyr::select, pathway, var1,
              var2,
              contrast, constant, padj, NES) %>%
          map(., ~.x %>%
                dplyr::mutate(direction = factor(ifelse(NES > 0, "up", "down")))) %>%
          map(droplevels) %>%
          map(., ~.x %>%
                dplyr::mutate(across(-c(padj, NES), as.character))))
  
  temp_list <- map(venn_list %>% names,
                   ~ setNames(venn_list[[.x]], paste0(names(venn_list[[.x]]), "_", .x))) %>%
    setNames(., names(venn_list))
  
  venn_df <- map(temp_list,
                 ~imap_dfr(., ~ .x %>%
                             bind_rows %>%
                             as_tibble %>%
                             dplyr::mutate(venn_name = .y) %>%
                             dplyr::select(venn_name, contrast, var1, 
                                           var2,
                                           constant, pathway, direction, NES, padj) %>%
                             pivot_longer(, cols = c(pathway),
                                          names_to = "metric",
                                          values_to = "pathway") %>%
                             dplyr::select(-metric) %>%
                             dplyr::mutate(across(-c(NES, padj), factor)) %>%
                             collect %>%
                             droplevels)) %>% 
    bind_rows(.)
  
  venn_list_plot <- venn_df %>% 
    split(., f = .$var1) %>%
    map(., ~.x %>%
          split(., f = .$direction) %>%
          map(., ~.x %>%
                split(., f = .$var2) %>%
                map(., ~.x %>%
                      dplyr::select(pathway) %>%
                      droplevels %>%
                      dplyr::mutate(across(everything(), as.character)) %>%
                      unlist) %>%
                Filter(function(k) length(k) > 0, .) #https://stackoverflow.com/questions/15292904/remove-empty-list-from-a-list-of-lists
          )) %>%
    map(., ~ .x %>%
          Filter(function(k) length(k) > 1, .))#have to remove contrasts in which var2 has only 1 factor
  
  venn_df <- venn_df %>%
    split(., f = .$contrast) %>%
    map(., ~.x %>% droplevels)
  
  allvars <- intersect((res_nociceptor_sig_list_renamed_tpm[["renamer"]] %>%
                          dplyr::select(cell_line, treatment, time) %>%
                          bind_rows %>%
                          distinct(.keep_all = TRUE) %>%
                          unlist %>%
                          as.character %>%
                          unique),
                       (imap_dfr(venn_df,
                                 ~ .x %>%
                                   as_tibble %>%
                                   dplyr::select(var1, var2) %>% 
                                   bind_rows %>% 
                                   distinct(.keep_all = TRUE)) %>%
                          unlist %>%
                          as.character %>%
                          unique))
  
  myColors <- viridisLite::turbo(n = length(allvars)) %>%
    setNames(., allvars[rank(allvars)])
  
  #for making the common legend:
  {
    venn_breaks_list <- F_venn_sde(venn_list_plot)
    
    temp_df <- venn_breaks_list[["legendBreaks"]] %>%
      map(., ~.x %>%
            unlist %>%
            max(.)) %>%
      unlist %>%
      which.max(.)
    
    venn_breaks_list <- list(
      list(venn_breaks_list[["myBreaks"]][[temp_df]], venn_breaks_list[["myBreaks"]][[temp_df]]) %>%
        setNames(., names(venn_df)),
      list(venn_breaks_list[["legendBreaks"]][[temp_df]], venn_breaks_list[["legendBreaks"]][[temp_df]]) %>%
        setNames(., names(venn_df))) %>%
      setNames(., c("myBreaks", "legendBreaks"))
    
    temp_df <- imap(venn_df,
                    ~ .x %>%
                      as_tibble %>%
                      dplyr::select(var1, var2) %>% 
                      pivot_longer(cols = -c(var1),
                                   names_to = "metric",
                                   values_to = "variable") %>%
                      dplyr::select(-metric) %>%
                      distinct(.keep_all = TRUE) %>%
                      collect %>%
                      droplevels)
    
    temp_df2 <- venn_breaks_list[["legendBreaks"]] %>%
      map(., ~.x %>%
            unlist %>%
            max(.) %>%
            tibble::enframe(., name = NULL, value = "max"))
    temp_df <- Map(c, temp_df, temp_df2)
    
    venn_legend <- map(temp_df %>% names,
                       ~list(
                         ((ggplot(data = temp_df[[.x]] %>% as.data.frame, aes(x = 1, y = 1, color = variable)) +
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
                         ((ggplot(data = data.frame(gene_count = c(0, temp_df[[.x]][["max"]])),
                                  aes(x = 1, y = 1, fill = gene_count)) +
                             geom_blank() +
                             geom_point() +
                             theme_void() +
                             scale_fill_viridis(discrete = FALSE, option = "mako",
                                                limits = c(0, temp_df[[.x]][["max"]]),
                                                name = "Significant pathways", guide = "colorbar")) %>%
                            g_legend())
                       ) %>%
                         setNames(., c("p1", "p2"))
    ) %>%
      setNames(., temp_df %>% names)
    }
  
  venn_breaks_list <- F_venn_sde(venn_list_plot)
  
  venn_fgsea_noci_direction <- list(venn_list_plot, venn_df, venn_breaks_list, venn_legend, allvars, myColors) %>%
    setNames(., c("venn_list_plot", "venn_df", "venn_breaks_list", "venn_legend", "allvars", "myColors"))
  
}
#venn_fgsea_noci_direction

#brca
{
  temp_df <- map(res_brca_fgsea_sig,
                 ~left_join(.x, res_brca_sig_list_renamed_tpm[["renamer"]] %>%
                              dplyr::select(label, cell_line, treatment, time) %>%
                              collect %>%
                              droplevels %>%
                              dplyr::mutate(across(.cols = everything(), as.factor))) %>% 
                   dplyr::select(renamed, label, cell_line, treatment, time, subset, pathway, padj, NES) %>% 
                   collect %>%
                   droplevels)
  
  #let's keep only "all_log2FC_tpm"
  venn_list <- temp_df[["all_log2FC_tpm"]] %>%
    bind_rows((dplyr::filter(., !grepl("_vs_", time)) %>%
                 dplyr::filter(., grepl("_vs_", treatment)) %>%
                 collect %>%
                 droplevels %>%
                 dplyr::mutate(contrast = "longitudinal",
                               var1 = time,
                               var2 = treatment,
                               constant = cell_line)),
              (dplyr::filter(., !grepl("_vs_", time)) %>%
                 dplyr::filter(., grepl("_vs_", treatment)) %>% #do you want a treatment-specific patch that has all timepoints for that treatment?
                 collect %>%
                 droplevels %>%
                 dplyr::mutate(contrast = "latitudinal",
                               var1 = treatment,
                               var2 = time,
                               constant = cell_line))) %>%  #timepoint-specific patch that has all treatments within that timepoint
    collect %>%
    droplevels %>%
    split(., f = .$contrast) %>%
    map(., droplevels) %>%
    map(., ~.x %>%
          split(., f = .$constant) %>%
          map(., ~.x %>%
                split(., f = .$var1) %>%
                map(as_tibble) %>%
                map(., dplyr::select, pathway, var1,
                    var2,
                    contrast, constant, padj, NES) %>%
                map(., ~.x %>%
                      dplyr::mutate(direction = factor(ifelse(NES > 0, "up", "down")))) %>%
                map(droplevels) %>%
                map(., ~.x %>%
                      dplyr::mutate(across(-c(padj, NES), as.character))))
    )
  
  temp_list <- map(venn_list %>% names,
                   ~ setNames(venn_list[[.x]], paste0(names(venn_list[[.x]]), "_", .x))) %>%
    setNames(., names(venn_list))
  
  venn_df <- map(temp_list,
                 ~imap_dfr(., ~ .x %>%
                             bind_rows %>%
                             as_tibble %>%
                             dplyr::mutate(venn_name = .y) %>%
                             dplyr::select(venn_name, contrast, var1, 
                                           var2,
                                           constant, pathway, direction, NES, padj) %>%
                             pivot_longer(, cols = c(pathway),
                                          names_to = "metric",
                                          values_to = "pathway") %>%
                             dplyr::select(-metric) %>%
                             dplyr::mutate(across(-c(NES, padj), factor)) %>%
                             collect %>%
                             droplevels)) %>% 
    bind_rows(.)
  # droplevels) %>%
  # split(., f = .$constant))
  # split(., f = .$venn_name))
  
  venn_list_plot <- venn_df %>% 
    split(., f = .$constant) %>%
    map(., ~.x %>%
          droplevels %>%
          split(., f = .$var1) %>%
          map(., ~.x %>%
                split(., f = .$direction) %>%
                map(., ~.x %>%
                      split(., f = .$var2) %>%
                      map(., ~.x %>%
                            dplyr::select(pathway) %>%
                            droplevels %>%
                            dplyr::mutate(across(everything(), as.character)) %>%
                            unlist) %>%
                      Filter(function(k) length(k) > 0, .) #https://stackoverflow.com/questions/15292904/remove-empty-list-from-a-list-of-lists
                ) 
              # )
          )) %>%
    map(., ~ .x %>%
          Filter(function(k) length(k) > 1, .))#have to remove contrasts in which var2 has only 1 factor
  
  #split it up: latitudinal in one list, longitudinal in other list
  
  venn_df <- venn_df %>% 
    split(., f = .$contrast) %>%
    map(., ~.x %>% droplevels) %>%
    map(., ~.x %>%
          split(., f = .$constant))
  
  venn_list_plot_lat <- map(venn_df[["latitudinal"]],
                            ~.x %>%
                              split(., f = .$direction) %>%
                              map(., ~.x %>%
                                    split(., f = .$var1) %>%
                                    map(., ~.x %>%
                                          split(., f = .$var2) %>%
                                          map(., ~.x %>%
                                                dplyr::select(pathway) %>%
                                                droplevels %>%
                                                dplyr::mutate(across(everything(), as.character)) %>%
                                                unlist) %>%
                                          Filter(function(k) length(k) > 0, .) #https://stackoverflow.com/questions/15292904/remove-empty-list-from-a-list-of-lists
                                    ) 
                              )
  ) %>%
    map(., ~ .x %>%
          Filter(function(k) length(k) > 1, .)) #have to remove contrasts in which var2 has only 1 factor
  
  venn_list_plot_long <- map(venn_df[["longitudinal"]],
                             ~.x %>%
                               split(., f = .$var1) %>%
                               map(., ~.x %>%
                                     split(., f = .$direction) %>%
                                     map(., ~.x %>%
                                           split(., f = .$var2) %>%
                                           map(., ~.x %>%
                                                 dplyr::select(pathway) %>%
                                                 droplevels %>%
                                                 dplyr::mutate(across(everything(), as.character)) %>%
                                                 unlist) %>%
                                           Filter(function(k) length(k) > 0, .) #https://stackoverflow.com/questions/15292904/remove-empty-list-from-a-list-of-lists
                                     ) )
  ) %>%
    map(., ~ .x %>%
          Filter(function(k) length(k) > 1, .)) #have to remove contrasts in which var2 has only 1 factor
  
  venn_breaks_list <- list(
    (map(venn_list_plot,
         ~map(., ~.x %>%
                purrr::flatten() %>%
                sapply(., length) %>%
                unlist %>%
                max(.) %>%
                signif(., digits = 2)) %>%
           map(., ~
                 c(seq(from = 0, to = .,
                       length.out = 101)))) %>%
       setNames(., names(venn_list_plot))),
    (map(venn_list_plot,
         ~map(., ~.x %>%
                purrr::flatten() %>%
                sapply(., length) %>%
                unlist %>%
                max(.) %>%
                signif(., digits = 2)) %>%
           map(., ~
                 c(seq(from = 0, to = max(.),
                       length.out = 7) %>% ceiling()
                 ))) %>%
       setNames(., names(venn_list_plot)))) %>%
    setNames(., c("myBreaks", "legendBreaks"))
  
  temp_df <- venn_breaks_list[["legendBreaks"]] %>%
    map(., ~.x %>%
          unlist %>%
          max(.) %>%
          tibble::enframe(., name = NULL, value = "max"))
  
  temp_df2 <- venn_df %>%
    purrr::flatten() %>%
    bind_rows() %>%
    split(., f = .$constant) %>%
    map(., ~.x %>%
          droplevels) %>%
    lapply(., get, x = "var1") %>%
    map(., ~.x %>%
          tibble::enframe(., name = NULL, value = "variable") %>%
          distinct())
  
  temp_df <- Map(c, temp_df, temp_df2)
  
  
  
  venn_df <- venn_df %>%
    purrr::flatten() %>%
    map(., ~.x %>%
          bind_rows() %>%
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
  
  venn_legend <- map(temp_df %>% names,
                     ~list(
                       ((ggplot(data = temp_df[[.x]] %>% as.data.frame, aes(x = 1, y = 1, color = variable)) +
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
                       ((ggplot(data = data.frame(gene_count = c(0, temp_df[[.x]][["max"]])),
                                aes(x = 1, y = 1, fill = gene_count)) +
                           geom_blank() +
                           geom_point() +
                           theme_void() +
                           scale_fill_viridis(discrete = FALSE, option = "mako",
                                              limits = c(0, temp_df[[.x]][["max"]]),
                                              name = "Significant pathways", guide = "colorbar")) %>%
                          g_legend())
                     ) %>%
                       setNames(., c("p1", "p2"))
  ) %>%
    setNames(., temp_df %>% names)
  
  venn_fgsea_brca_direction <- list(venn_list_plot_lat, venn_list_plot_long, venn_df, venn_breaks_list, venn_legend, allvars, myColors) %>%
    setNames(., c("venn_list_plot_lat", "venn_list_plot_long", "venn_df", "venn_breaks_list", "venn_legend", "allvars", "myColors"))
  
}
#venn_fgsea_brca_direction




write_rds(venn_fgsea_noci_direction, paste0(dir, "/data/venn_fgsea_noci_direction.rds"), compress = "gz")
write_rds(venn_fgsea_brca_direction, paste0(dir, "/data/venn_fgsea_brca_direction.rds"), compress = "gz")


# Make Venn figures -------------------------------------------------------

venn_list <- map(venn_fgsea_noci_direction[["venn_df"]],
                 ~ .x %>%
                   as_tibble %>%
                   dplyr::select(var1) %>% 
                   distinct(.keep_all = TRUE)) %>%
  map(., ~.x %>%
        unlist(.) %>%
        as.character %>%
        map(., 
            ~((patch_venn2(venn_fgsea_noci_direction[["venn_list_plot"]][[.x]],
                           venn_fgsea_noci_direction[["myColors"]],
                           venn_fgsea_noci_direction[["venn_breaks_list"]][["myBreaks"]][[.x]]) +
                 theme(legend.position = "none") +
                 # plot_layout(ncol = 2)))) %>%
                 plot_layout(ncol = 2))) +
              plot_annotation(tag_levels = list(.x))) %>%
        # plot_layout(ncol = 2))) + ggtitle(.x)) %>%
        setNames(., unlist(.x)))


map(venn_list %>% names,
    ~lapply(venn_list[[.x]] %>% names, function(y) {
      venn_list[[.x]][[y]] %>%
        wrap_elements(.) + ggtitle(y)
    }) %>%
      purrr::reduce(., `/`) + plot_layout(nrow = 3) | ((wrap_elements(venn_fgsea_noci_direction[["venn_legend"]][[.x]][["p2"]], ignore_tag = TRUE) / wrap_elements(venn_fgsea_noci_direction[["venn_legend"]][[.x]][["p1"]], ignore_tag = TRUE))) +
      plot_layout(widths = c(3,1))) %>%
  setNames(., names(venn_list)) %>%
  map2(., names(.),
       ~ggsave(paste0(dir, "/figures/", Sys.Date(), "-venn_noci_fgsea_treat_dir-", .y, ".png"),
               .x,
               width = 8, height = 8, units = "in"))


# Venn figure for brca ----------------------------------------------------


venn_list <- map(venn_fgsea_brca_direction[["venn_list_plot_lat"]] %>% names,
                 ~map2(., names(venn_fgsea_brca_direction[["venn_list_plot_lat"]][[.x]]),
                       ~((patch_venn2(venn_fgsea_brca_direction[["venn_list_plot_lat"]][[.x]][[.y]],
                                      venn_fgsea_brca_direction[["myColors"]],
                                      venn_fgsea_brca_direction[["venn_breaks_list"]][["myBreaks"]][[.x]] %>% unlist) +
                            theme(legend.position = "none") +
                            plot_layout(ncol = 1))) +
                         plot_annotation(tag_levels = list(.y))
                 ) %>%
                   setNames(., names(venn_fgsea_brca_direction[["venn_list_plot_lat"]][[.x]]))
) %>%
  setNames(., venn_fgsea_brca_direction[["venn_list_plot_lat"]] %>% names)

map(venn_list %>% names,
    ~((wrap_elements(venn_list[[.x]][["down"]]) | wrap_elements(venn_list[[.x]][["up"]])) |
        (wrap_elements(venn_fgsea_brca_direction[["venn_legend"]][[.x]][["p2"]], ignore_tag = TRUE) / wrap_elements(venn_fgsea_brca_direction[["venn_legend"]][[.x]][["p1"]], ignore_tag = TRUE))) +
      plot_layout(widths = c(2,2,1)) +
      plot_annotation(title = paste0("Direction of significant pathways within ", paste0(.x)))) %>%
  setNames(., names(venn_list)) %>%
  map2(., names(.),
       ~ggsave(paste0(dir, "/figures/", Sys.Date(), "-venn_brca_fgsea_treat_dir-", .y, ".png"),
               .x,
               width = 10, height = 7, units = "in"))

venn_list <- map(venn_fgsea_brca_direction[["venn_list_plot_long"]] %>% names,
                 ~map2(., names(venn_fgsea_brca_direction[["venn_list_plot_long"]][[.x]]),
                       ~((patch_venn2(venn_fgsea_brca_direction[["venn_list_plot_long"]][[.x]][[.y]],
                                      venn_fgsea_brca_direction[["myColors"]],
                                      venn_fgsea_brca_direction[["venn_breaks_list"]][["myBreaks"]][[.x]] %>% unlist) +
                            theme(legend.position = "none") +
                            plot_layout(ncol = 2))) +
                         plot_annotation(tag_levels = list(.y))
                 ) %>%
                   setNames(., names(venn_fgsea_brca_direction[["venn_list_plot_long"]][[.x]]))
) %>%
  setNames(., venn_fgsea_brca_direction[["venn_list_plot_long"]] %>% names)

map(venn_list %>% names,
    ~((wrap_elements(venn_list[[.x]][["1h"]]) / wrap_elements(venn_list[[.x]][["6h"]]) / wrap_elements(venn_list[[.x]][["24h"]])) |
        (wrap_elements(venn_fgsea_brca_direction[["venn_legend"]][[.x]][["p2"]], ignore_tag = TRUE) / wrap_elements(venn_fgsea_brca_direction[["venn_legend"]][[.x]][["p1"]], ignore_tag = TRUE))) +
      plot_layout(widths = c(2,1)) +
      plot_annotation(title = paste0("Direction of significant pathways within ", paste0(.x)))) %>%
  setNames(., names(venn_list)) %>%
  map2(., names(.),
       ~ggsave(paste0(dir, "/figures/", Sys.Date(), "-venn_brca_fgsea_time_dir-", .y, ".png"),
               .x,
               width = 10, height = 7, units = "in"))


