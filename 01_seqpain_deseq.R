## Sharon Grim
## 2023

## 01_seqpain_deseq.R

##This script takes the "counts" tables of genes in samples,
##and runs DESeq analyses on them.
##All scripts assume we are working in a local directory called "seqpain"
##and that data files are within "/seqpain/data/"
##For ease, I preprocessed AnnotationHub reference data: gene_tx_EnsDb.rds


library(tidyverse)
library(DESeq2)
library(BiocParallel)

##if you have a modern computer, you can use multicore processing.
cluster <- multidplyr::new_cluster(4)
multidplyr::cluster_library(cluster, "dplyr")
bpparam_multi <- MulticoreParam(timeout = 100,
                                workers = 4,
                                stop.on.error = TRUE,
                                RNGseed = 48105,
                                progressbar = TRUE)
register(bpparam_multi)

dir <- getwd()

##but if you don't, please check the functions and comment out where indicated.

##helper function:
coalesce2 <- function(x, y, sep = ".") ifelse(x == y, coalesce(x, y, sep = sep), paste0(x, "_vs_", y))



# Custom functions --------------------------------------------------------

#This performs DESeq analyses (multicore enabled)
F_deseq <- function(counts_df, mdata, cell, deseq_design){
  #d
  deseq_design <- rlang::parse_expr(deseq_design)
  cell <- rlang::parse_expr(cell)
  #m
  metadata <- mdata %>%
    dplyr::filter(cell_type %in% c({cell})) %>%
    droplevels(.) %>%
    mutate(treatment = factor(treatment, levels = c("DMSO", "Kw", "Taxol", "KwTaxol")),
           cell_line = factor(cell_line),
           cell_type = factor(cell_type),
           SeqRun = factor(SeqRun)) %>%
    mutate(timeindex = as.character(time) %>%
             fct_inseq() %>%
             as.integer(),
           time = factor(time),
           timepoint = factor(timepoint))
  #data
  matrix <- counts_df %>%
    dplyr::filter(sample_name %in% levels(metadata$sample_name)) %>%
    ungroup() %>%
    dplyr::select(gene_symbol, sample_name, counts) %>%
    droplevels(.) %>%
    mutate(across(!contains("count"), factor)) %>%
    group_by(sample_name, gene_symbol) %>%
    dplyr::summarise(Reads = sum(counts, na.rm = TRUE), .groups = "keep") %>%
    drop_na() %>%
    pivot_wider(., id_cols = "gene_symbol",
                names_from = "sample_name",
                values_fill = 0,
                values_from = "Reads") %>%
    column_to_rownames(., var = "gene_symbol")
  
  dds <- DESeqDataSetFromMatrix(countData = matrix,
                                tidy = FALSE,
                                colData = metadata,
                                design = eval(deseq_design)) %>%
    DESeq(., 
          # parallel = FALSE)
          parallel = TRUE, BPPARAM = bpparam_multi) #comment this out if you do not have multicore processing enabled
  
  mod_mat <- model.matrix(eval(deseq_design),
                          data = metadata) %>%
    as.data.frame(.) %>%
    dplyr::rename_with(., ~ gsub(":", "_vs_", .x, fixed = TRUE),  contains(":")) %>%
    as.matrix(.)
  
  list(metadata, matrix, dds, mod_mat) %>%
    setNames(., c("metadata", "matrix", "dds", "model_matrix"))
  
}

#Generate a DESeq results list (multicore enabled)
F_extract_results_multi <- function(dlist, set.alpha){
  de <- dlist[["dds"]]
  metadata <- dlist[["metadata"]]
  
  #when there is a 0 intercept, then the labeling of contrast is "suchandsuchA" instead of "suchandsuch_A"
  rn <- resultsNames(de) %>%
    setdiff(c("0", "Intercept")) %>%
    as.data.frame() %>%
    setNames(., "label")
  
  #for the extra contrasts:
  #remove the spurious comparisons (e.g. treatmentA_celllineA vs treatmentB_celllineB)
  drop <- as.character(design(de)) %>%
    stringr::str_replace_all(., "[\\s~]", "") %>%
    strsplit(., "[\\*+]", perl = TRUE) %>%
    unlist(.) %>%
    setdiff(c("0", "Intercept"))
  
  for(i in 1:length(drop)){
    if(length(levels(metadata[[drop[i]]])) <= 2){ #if there are not enough levels in this metadata factor to make further combinations:
      temp_df <- paste0(drop[i], "_", levels(metadata[[drop[i]]])[-1], "_vs_", levels(metadata[[drop[i]]])[1]) %>%
        as.data.frame() %>%
        setNames(., "label")
      if(exists("rn")){
        rn <- bind_rows(rn, temp_df)}
      if(!exists("rn")){
        rn <- temp_df}
    } else {
      temp_df <- paste0(drop[i], "_", levels(metadata[[drop[i]]])[-1], "_vs_", levels(metadata[[drop[i]]])[1]) %>%
        rev() %>%
        combn(., 2) %>%
        t() %>%
        as.data.frame() %>%
        mutate(label = paste0(str_replace(.[,1], "(_vs_).+$", ""), "_vs_", str_replace(.[,2], "(_vs_).+$", ""))) %>%
        mutate(label = str_replace(label, paste0("_vs_", drop[i], "_"), "_vs_"))
      if(exists("rn2")){
        rn2 <- bind_rows(rn2, temp_df)}
      if(!exists("rn2")){
        rn2 <- temp_df}
    } 
  }
  
  #this is for multifactor contrasts
  if(length(drop) - 1 > 0){
    for(i in 1:(length(drop)-1)){
      if(length(levels(metadata[[drop[i+1]]])) > 2){
        temp_df <- paste0(drop[i+1], levels(metadata[[drop[i+1]]])[-1]) %>%
          rev() %>%
          combn(., 2) %>%
          t() %>%
          as.data.frame() %>%
          dplyr::slice(rep(1:n(), each = length(levels(metadata[[drop[i]]])[-1]))) %>%
          mutate(across(everything(), ~paste0(drop[i], levels(metadata[[drop[i]]])[-1], ".", .x)))
        temp_df <- temp_df %>%
          mutate(label = paste0(str_replace(.[,1], "(_vs_).+$", ""), "_vs_", str_replace(.[,2], "(_vs_).+$", ""))) %>%
          mutate(label = str_replace(label, "(_vs_).+\\.", "_vs_"))
        if(exists("rn2")){
          rn2 <- bind_rows(rn2, temp_df)}
        if(!exists("rn2")){
          rn2 <- temp_df}
      }
    }
  }
  
  #remove the spurious comparisons (e.g. treatmentA_celllineA vs treatmentB_celllineB)
  rename_labels <- map(drop %>% set_names(),
                       ~dplyr::filter(rn2 %>%
                                        distinct(label) %>%
                                        mutate(type = label), grepl(.x, label)) %>%
                         mutate(pair = gsub(paste0(.x, "_"), "", label)) %>%
                         separate(pair, into = c("first", "second"), sep = "_vs_", extra = "merge", remove = TRUE) %>%
                         mutate(first = factor(first),
                                second = factor(second)) %>%
                         as.data.frame() %>%
                         left_join(., metadata %>%
                                     dplyr::select(all_of(.x), cell_type, cell_line, treatment, time, SeqRun) %>%
                                     droplevels(.) %>%
                                     mutate(time = paste0(time, "h")),
                                   by = c("first" = .x)) %>%
                         left_join(., metadata %>%
                                     dplyr::select(all_of(.x), cell_type, cell_line, treatment, time, SeqRun) %>%
                                     droplevels(.) %>%
                                     mutate(time = paste0(time, "h")),
                                   suffix = c(".first", ".second"),
                                   by = c("second" = .x))) %>%
    bind_rows(., .id = "type") %>%
    mutate(cell_type = across(matches("cell_type")) %>%
             purrr::reduce(coalesce2),
           cell_line = across(matches("cell_line")) %>%
             purrr::reduce(coalesce2),
           treatment = across(contains("treatment")) %>%
             purrr::reduce(coalesce2),
           time = across(matches("^time")) %>%
             purrr::reduce(coalesce2),
           SeqRun = across(matches("SeqRun")) %>%
             purrr::reduce(coalesce2)) %>%
    mutate(across(c(cell_type, cell_line, treatment, time, SeqRun), factor))
  
  keep <- combn(c("cell_type", "cell_line", "treatment", "time", "SeqRun"), 4) %>%
    as.data.frame() 
  #tidy the combo contrasts dataframe
  rn2 <- map(keep %>% set_names(),
             ~dplyr::filter(rename_labels, !Reduce(`|`, across(all_of(.x), ~ grepl("_vs_", .x)))) %>%
               as.data.frame(.)) %>%
    bind_rows(.) %>%
    dplyr::select(label) %>%
    droplevels(.) %>%
    left_join(rn2) %>%
    relocate(-label) %>%
    distinct(label, .keep_all = TRUE)
  
  rn <- rn %>%
    distinct()
  
  #for single factor contrasts
  if(exists("rn")){
    for(i in 1:nrow(rn)){
      res1 <- results(de,
                      contrast = list(rn[i, 1]),
                      # parallel=FALSE,
                      parallel=TRUE, BPPARAM=bpparam_multi, #multicore processing
                      alpha = set.alpha) %>%
        as.data.frame() %>%
        rownames_to_column("gene_symbol") %>%
        left_join(., (lfcShrink(de,
                                contrast = list(rn[i, 1]),
                                type = "ashr",
                                # parallel=FALSE,
                                parallel=TRUE, BPPARAM=bpparam_multi, #multicore processing
                                quiet = TRUE) %>%
                        as.data.frame() %>%
                        rownames_to_column("gene_symbol") %>%
                        dplyr::rename(log2FoldChange_MMSE = log2FoldChange) %>%
                        dplyr::select(gene_symbol, log2FoldChange_MMSE)),
                  by = "gene_symbol") %>%
        mutate(contrast = as.character(rn[i, 1])) %>%
        relocate(contrast, .before = "gene_symbol")
      if(exists("res")){
        res <- bind_rows(res, res1)}
      if(!exists("res")){
        res <- res1}
    }
  }
  
  #for multifactor contrasts
  if(exists("rn2")){
    for(i in 1:nrow(rn2)){
      res1 <- results(de,
                      contrast = list(rn2[i, 1], rn2[i, 2]),
                      # parallel=FALSE,
                      parallel=TRUE, BPPARAM=bpparam_multi, #multicore processing
                      alpha = set.alpha) %>%
        as.data.frame() %>%
        rownames_to_column("gene_symbol") %>%
        left_join(., lfcShrink(de,
                               contrast = list(rn2[i, 1], rn2[i, 2]),
                               # parallel=FALSE,
                               parallel=TRUE, BPPARAM=bpparam_multi, #multicore processing
                               type = "ashr",
                               quiet = TRUE) %>%
                    as.data.frame() %>%
                    rownames_to_column("gene_symbol") %>%
                    dplyr::rename(log2FoldChange_MMSE = log2FoldChange) %>%
                    dplyr::select(gene_symbol, log2FoldChange_MMSE),
                  by = "gene_symbol") %>%
        mutate(contrast = rn2[i, "label"]) %>%
        relocate(contrast, .before = "gene_symbol")
      if(exists("res")){
        res <- bind_rows(res, res1)}
      if(!exists("res")){
        res <- res1}
    }
  }
  bind_rows(res) %>%
    mutate(across(c(contrast, gene_symbol), factor))
}

#This pulls out all results from DEseq analyses with p < set.alpha (multicore enabled)
F_keep_sde <- function(resultslist, headers, set.alpha){
  headers <- intersect(headers, names(resultslist))
  
  df <- resultslist %>%
    dplyr::select(all_of(headers)) %>%
    dplyr::arrange(by = desc(-padj)) %>%
    dplyr::mutate(contrast = factor(contrast)) %>%
    dplyr::mutate(padj_filtered = replace(padj, which(padj > set.alpha), NA)) %>%
    dplyr::slice(which(!is.na(padj_filtered))) %>%
    dplyr::select(-padj_filtered) %>%
    split(f = .$contrast, drop = TRUE) %>%
    lapply(., function(x) {
      x %>% dplyr::select(!c(contrast))})
  
  namevar <- deparse(substitute(resultslist))
  assign(paste0(namevar, "_sig_list"), df, envir = .GlobalEnv)
  
  df2 <- df %>%
    imap_dfr(., ~
               .x %>%
               mutate(contrast = .y) %>%
               pivot_longer(cols = -c(contrast, gene_symbol, baseMean, padj),
                            names_to = "metric",
                            values_to = "log2FoldChange_MMSE")) %>%
    dplyr::select(-metric) %>%
    droplevels(.) %>%
    mutate(across(c(gene_symbol, contrast), factor)) %>%
    multidplyr::partition(cluster) %>% #comment out if not using multicore
    group_by(contrast) %>%
    dplyr::summarise(num_sig_genes = length(gene_symbol)) %>%
    collect %>%
    droplevels %>%
    mutate(num_all_genes = length(unique(resultslist$gene_symbol)),
           proportion = num_sig_genes / num_all_genes * 100)
  
  namevar <- deparse(substitute(resultslist))
  assign(paste0("summary_", namevar, "_sig"), df2, envir = .GlobalEnv)
  rm(df)
  rm(df2)
}

#This calculates TPM for each gene, and is used in F_rename_sde
F_calculateTPM <- function(cell, mdata, counts_df){
  cell <- rlang::parse_expr(cell)
  metadata <- mdata %>%
    dplyr::filter(cell_type %in% c({cell})) %>%
    droplevels(.)
  #data
  tpm <- full_join(
    (counts_df %>%
       dplyr::filter(sample_name %in% levels(metadata$sample_name)) %>%
       ungroup() %>%
       dplyr::select(sample_name, counts) %>%
       droplevels(.) %>%
       mutate(across(!contains("count"), factor)) %>%
       group_by(sample_name) %>%
       dplyr::summarise(sumA = sum(counts, na.rm = TRUE), .groups = "keep") %>%
       collect),
    (counts_df %>%
       dplyr::filter(sample_name %in% levels(metadata$sample_name)) %>%
       ungroup() %>%
       dplyr::select(c(gene_symbol, sample_name, counts)) %>%
       droplevels(.) %>%
       mutate(across(!contains("count"), factor)) %>%
       group_by(gene_symbol, sample_name) %>%
       dplyr::summarise(A = sum(counts, na.rm = TRUE), .groups = "keep") %>%
       collect)) %>%
    collect %>%
    droplevels %>%
    # multidplyr::partition(cluster) %>% #comment out if not using multicore
    # group_by(gene_symbol, sample_name) %>%
    mutate(TPM = A / sumA * 1000000) %>%
    # collect %>%
    ungroup() %>%
    dplyr::select(-c(A, sumA)) %>%
    droplevels() %>%
    # distinct(sample_name, gene_symbol, .keep_all = TRUE) %>%
    mutate(across(!contains("TPM"), factor)) %>%
    mutate(sample_name = factor(sample_name, levels = levels(rnaseq_metadata$sample_name)))
}

#This renames the "groups" to something more legible, pulls out two abundance metrics for each gene (avgTPM and DESeq-normalized abundance)
#and creates a nice list of those results (multicore enabled)
F_rename_sde <- function(dlist, rlist){
  de <- dlist[["dds"]]
  counts <- dlist[["matrix"]]  %>%
    rownames_to_column(var = "gene_symbol") %>%
    pivot_longer(., cols = !c(gene_symbol),
                 names_to = "sample_name",
                 values_to = "counts")
  metadata <- dlist[["metadata"]]
  namevar <- deparse(substitute(rlist))
  cell <- deparse(substitute(dlist)) %>%
    gsub("(deseq_)|(res_)", "", .)
  
  temp_df <- NULL
  rn2 <- NULL
  rn <- resultsNames(de) %>%
    setdiff(c("0", "Intercept")) %>%
    as.data.frame() %>%
    setNames(., "label")
  
  #for the extra contrasts:
  #remove the spurious comparisons (e.g. treatmentA_celllineA vs treatmentB_celllineB)
  drop <- as.character(design(de)) %>%
    stringr::str_replace_all(., "[\\s~]", "") %>%
    strsplit(., "[\\*+]", perl = TRUE) %>%
    unlist(.) %>%
    setdiff(c("0", "Intercept"))
  
  tpm_df <- F_calculateTPM(cell, metadata, counts)
  
  tpm_avg_by_contrast <- map(drop %>% set_names(),
                             ~dplyr::select(metadata, 
                                            c(all_of(.x), sample_name)) %>%
                               droplevels(.) %>%
                               as.data.frame() %>%
                               left_join(., (tpm_df %>%
                                               dplyr::filter(sample_name %in% levels(metadata$sample_name)) %>%
                                               ungroup() %>%
                                               droplevels(.) %>%
                                               mutate(across(!contains("TPM"), factor))),
                                         by = "sample_name") %>%
                               droplevels(.) %>%
                               mutate(across(c(all_of(.x), gene_symbol), factor),
                                      sample_name = as.character(sample_name)) %>%
                               # multidplyr::partition(cluster) %>% #comment out if not multicore processing
                               group_by(across(c(all_of(.x), gene_symbol))) %>%
                               # group_by(across(where(is.factor))) %>%
                               dplyr::summarise(avgTPM = mean(TPM, na.rm = TRUE),
                                                sdTPM = sd(TPM, na.rm = TRUE),
                                                log2avgTPM = log2(avgTPM+1), .groups = "keep") %>%
                                                # log2avgTPM = log2(avgTPM+1)) %>%
                               # collect %>%
                               droplevels(.)) %>%
    bind_rows(., .id = "type") %>%
  dplyr::mutate(type = factor(type))
  
  for(i in 1:length(drop)){
    if(length(levels(metadata[[drop[i]]])) <= 2){ #if there are not enough levels in this metadata factor to make further combinations:
      temp_df <- paste0(drop[i], "_", levels(metadata[[drop[i]]])[-1], "_vs_", levels(metadata[[drop[i]]])[1]) %>%
        as.data.frame() %>%
        setNames(., "label")
      if(exists("rn")){
        rn <- bind_rows(rn, temp_df)}
      if(!exists("rn")){
        rn <- temp_df}
    } else {
      temp_df <- paste0(drop[i], "_", levels(metadata[[drop[i]]])[-1], "_vs_", levels(metadata[[drop[i]]])[1]) %>%
        rev() %>%
        combn(., 2) %>%
        t() %>%
        as.data.frame() %>%
        mutate(label = paste0(str_replace(.[,1], "(_vs_).+$", ""), "_vs_", str_replace(.[,2], "(_vs_).+$", ""))) %>%
        mutate(label = str_replace(label, paste0("_vs_", drop[i], "_"), "_vs_"))
      if(exists("rn2")){
        rn2 <- bind_rows(rn2, temp_df)}
      if(!exists("rn2")){
        rn2 <- temp_df}
    } 
  }
  
  #this is for multifactor contrasts
  if(length(drop) - 1 > 0){
    for(i in 1:(length(drop)-1)){
      if(length(levels(metadata[[drop[i+1]]])) > 2){
        temp_df <- paste0(drop[i+1], levels(metadata[[drop[i+1]]])[-1]) %>%
          rev() %>%
          combn(., 2) %>%
          t() %>%
          as.data.frame() %>%
          dplyr::slice(rep(1:n(), each = length(levels(metadata[[drop[i]]])[-1]))) %>%
          mutate(across(everything(), ~paste0(drop[i], levels(metadata[[drop[i]]])[-1], ".", .x)))
        temp_df <- temp_df %>%
          mutate(label = paste0(str_replace(.[,1], "(_vs_).+$", ""), "_vs_", str_replace(.[,2], "(_vs_).+$", ""))) %>%
          mutate(label = str_replace(label, "(_vs_).+\\.", "_vs_"))
        if(exists("rn2")){
          rn2 <- bind_rows(rn2, temp_df)}
        if(!exists("rn2")){
          rn2 <- temp_df}
      }
    }
  }
  rn <- rn %>%
    distinct()
  
  #remove the spurious comparisons (e.g. treatmentA_celllineA vs treatmentB_celllineB)
  rename_labels <- map(drop %>% set_names(),
                       ~dplyr::filter(rn2 %>%
                                        distinct(label) %>%
                                        bind_rows(., rn) %>%
                                        mutate(type = label), grepl(.x, label)) %>%
                         mutate(pair = gsub(paste0(.x, "_"), "", label)) %>%
                         separate(pair, into = c("first", "second"), sep = "_vs_", extra = "merge", remove = TRUE) %>%
                         mutate(first = factor(first),
                                second = factor(second)) %>%
                         as.data.frame() %>%
                         left_join(., metadata %>%
                                     ungroup() %>%
                                     dplyr::select(c(all_of(.x), cell_type, cell_line, treatment, time, SeqRun)) %>%
                                     droplevels(.) %>%
                                     mutate(time = paste0(time, "h"),
                                            FC = SeqRun) %>%
                                     mutate(across(everything(), factor)) %>%
                                     mutate(time = factor(time, levels = c("0h", "1h", "6h", "24h")),
                                            treatment = factor(treatment, levels = c("DMSO", "Kw", "Taxol", "KwTaxol"))) %>%
                                     droplevels,
                                   by = c("first" = .x)) %>%
                         left_join(., metadata %>%
                                     ungroup() %>%
                                     dplyr::select(c(all_of(.x), cell_type, cell_line, treatment, time, SeqRun)) %>%
                                     droplevels(.) %>%
                                     mutate(time = paste0(time, "h"),
                                            FC = SeqRun) %>%
                                     mutate(across(everything(), factor)) %>%
                                     mutate(time = factor(time, levels = c("0h", "1h", "6h", "24h")),
                                            treatment = factor(treatment, levels = c("DMSO", "Kw", "Taxol", "KwTaxol"))) %>%
                                     droplevels,
                                   suffix = c(".first", ".second"),
                                   by = c("second" = .x))) %>%
    bind_rows(., .id = "type") %>%
    mutate(cell_type = across(matches("cell_type")) %>%
             purrr::reduce(coalesce2),
           cell_line = across(matches("cell_line")) %>%
             purrr::reduce(coalesce2),
           treatment = across(contains("treatment")) %>%
             purrr::reduce(coalesce2),
           time = across(matches("^time")) %>%
             purrr::reduce(coalesce2),
           # SeqRun = across(matches("SeqRun")) %>%
           FC = across(matches("FC")) %>%
             purrr::reduce(coalesce2)) %>%
    mutate(across(c(cell_type, cell_line, treatment, time, FC), factor)) %>%
    distinct(.)
  
  # 
  #tidy the combo contrasts dataframe
  renamer <- bind_rows(
    (map(combn(c("treatment", "time"), 1) %>% #remove any rows where the contrast is between different times and different treatments
           as.data.frame() %>% set_names(),
         ~dplyr::filter(rename_labels %>%
                          dplyr::slice(which(!grepl("SeqRun", type))) %>%
                          dplyr::slice(which(!grepl("_vs_", cell_line))) %>% #remove rows where the contrast is between different cell lines
                          distinct(), !Reduce(`|`, across(all_of(.x), ~ grepl("_vs_", .x)))) %>%
           as.data.frame(.)) %>%
       bind_rows() %>%
       mutate(renamed = paste(treatment, cell_line, time, sep = ":")) %>%
       distinct(renamed, .keep_all = TRUE)),
    ( map(combn(c("treatment", "time", "FC"), 3) %>% #contrast is between different cell lines at exact same treatment, time, and SeqRun
            as.data.frame() %>% set_names(),
          ~dplyr::filter(rename_labels %>%
                           dplyr::slice(which(!grepl("SeqRun", type))) %>%
                           dplyr::slice(which(grepl("_vs_", cell_line))) %>% #keep rows where the contrast is between different cell lines
                           distinct(), !Reduce(`|`, across(all_of(.x), ~ grepl("_vs_", .x)))) %>%
            as.data.frame(.)) %>%
        bind_rows() %>%
        mutate(renamed = paste(treatment, cell_line, time, sep = ":")) %>%
        distinct(renamed, .keep_all = TRUE)),
    (map(
      combn(c("cell_type", "cell_line", "FC"), 3) %>% #remove any rows where the contrast is between 3 of these 4 
        as.data.frame() %>% set_names(),
      ~dplyr::filter(rename_labels %>%
                       dplyr::slice(which(!grepl("SeqRun", type))) %>%
                       dplyr::slice(which(grepl("DMSO", treatment.second))), !Reduce(`|`, across(.x, ~ grepl("_vs_", .x)))) %>%
        as.data.frame(.)) %>%
       bind_rows(.) %>%
       distinct() %>%
       mutate(renamed = paste(treatment, cell_line, time, sep = ":"))),
    (map(
      combn(c("cell_type", "cell_line", "treatment", "time"), 4) %>% #remove any rows where the contrast is between 3 of these 4 
        as.data.frame() %>% set_names(),
      ~dplyr::filter(rename_labels %>%
                       dplyr::slice(which(grepl("SeqRun", type))), !Reduce(`|`, across(.x, ~ grepl("_vs_", .x)))) %>%
        as.data.frame(.)) %>%
       bind_rows(.) %>%
       distinct() %>%
       mutate(renamed = paste(FC, treatment, cell_line, time, sep = ":"))),
    (rename_labels %>%
       dplyr::filter(label %in% rn$label) %>%
       ungroup() %>%
       distinct() %>%
       mutate(renamed = paste(treatment, cell_line, time, sep = ":")))) %>%
    dplyr::select(label, renamed, type, first, second, cell_type, cell_line, treatment, time, FC) %>%
    dplyr::mutate(renamed = ifelse(type == "SeqRun", label, renamed)) %>%
    distinct(renamed, .keep_all = TRUE) 
  
  rn <- rn %>%
    distinct() %>%
    left_join(rename_labels) %>%
    tidyr::drop_na(label) %>%
    mutate(renamed = ifelse(type == "SeqRun", label, paste(treatment, cell_line, time, sep = ":"))) %>%
    dplyr::select(label, renamed) %>%
    relocate(-c(label, renamed)) %>%
    distinct(label, .keep_all = TRUE)
  
  rn2 <- renamer %>%
    dplyr::select(label, renamed, first, second, cell_type, cell_line, treatment, time, FC) %>%
    left_join(rn2) %>%
    tidyr::drop_na() %>%
    anti_join(rn, by = c("label", "renamed")) %>%
    dplyr::select(label, renamed) %>%
    relocate(-c(label, renamed)) %>%
    distinct(label, renamed, .keep_all = TRUE)
  
  renamer <- bind_rows(rn, rn2) %>%
    left_join(renamer) %>%
    droplevels
  
  #Fix the "renamer" df so that it makes more sense.
  renamer <- renamer %>%
    dplyr::mutate(type = factor(type, levels = sort(unique(type))),
                  time = factor(time, levels = c("0h", "1h_vs_0h",
                                                 "1h", "1h_vs_6h", "1h_vs_24h",
                                                 "6h_vs_0h", "6h_vs_1h", "6h", "6h_vs_24h",
                                                 "24h_vs_0h", "24h_vs_1h", "24h_vs_6h", "24h")),
                  treatment = factor(treatment, levels = c("DMSO", "Kw", "KwTaxol", "Taxol", 
                                                           "Kw_vs_DMSO", "KwTaxol_vs_DMSO", "Taxol_vs_DMSO", 
                                                           "KwTaxol_vs_Taxol", "KwTaxol_vs_Kw", "Taxol_vs_Kw")),
                  cell_line = factor(cell_line, levels = c("iPSCs", "Hs578T", "MDAMB231", "MDAMB231_vs_Hs578T", "MDAMB231_vs_MCF7", 
                                                           "MCF7", "MCF7_vs_Hs578T"))) %>%
    arrange(., by = type, cell_line, treatment, time) %>%
    dplyr::mutate(renamed = fct_relevel(renamed, .$renamed[order(.$renamed, .$cell_line, .$treatment, .$time)])) %>% #for brca
    dplyr::mutate(renamed = factor(renamed, levels = unique(renamed))) %>%
    collect %>%
    droplevels
  
  res_sig_renamed_avgtpm <- renamer %>%
    dplyr::select(label, renamed) %>%
    distinct(.) %>%
    right_join(., (rlist %>%
                     imap_dfr(., ~
                                .x %>%
                                mutate(contrast = .y) %>%
                                pivot_longer(cols = -c(contrast, gene_symbol, baseMean, padj),
                                             names_to = "metric",
                                             values_to = "log2FoldChange_MMSE")) %>%
                     dplyr::select(-metric) %>%
                     droplevels(.) %>%
                     mutate(across(c(gene_symbol, contrast), factor)) %>%
                     bind_rows(.) %>%
                     ungroup(.)),
               by = c("label" = "contrast")) %>%
    left_join(., rename_labels %>%
                dplyr::select(label, type, first, second) %>%
                distinct(.)) %>%
    dplyr::mutate(across(!c(baseMean, padj, log2FoldChange_MMSE), factor))
  
  
  res_list <- list(tpm_avg_by_contrast, res_sig_renamed_avgtpm, renamer) %>%
    setNames(., c("tpm_avg_by_contrast", "res_sig_renamed_avgtpm", "renamer"))
  assign(paste0(namevar, "_renamed_tpm"), res_list, envir = .GlobalEnv)
}

#This creates a long df of the significant results from DESeq. If you want it.
F_extract_sde <- function(dlist, rlist){
  res_list <- F_rename_labels(dlist, rlist)
  
  drop <- as.character(design(dlist[["dds"]])) %>%
    stringr::str_replace_all(., "[\\s~]", "") %>%
    strsplit(., "[\\*+]", perl = TRUE) %>%
    unlist(.)
  
  df <- map(drop %>% set_names(),
            ~dplyr::filter(res_list[["res_sig_renamed_avgtpm"]] %>% 
                             dplyr::select(gene_symbol, label, renamed, log2FoldChange_MMSE, type, first, second),
                           grepl(.x, type)) %>%
              pivot_longer(., cols = c(first, second),
                           names_to = "metric",
                           values_to = "group") %>%
              ungroup() %>%
              dplyr::select(-c(metric, type)) %>%
              droplevels(.) %>%
              mutate(across(c(gene_symbol, label, renamed, group), factor)) %>%
              distinct(.) %>%
              left_join(., res_list[["tpm_avg_by_contrast"]] %>%
                          ungroup(.) %>%
                          dplyr::select(c(all_of(.x), gene_symbol, log2avgTPM)) %>%
                          droplevels() %>%
                          pivot_longer(., cols = c(all_of(.x)),
                                       names_to = "metric",
                                       values_to = "group") %>%
                          dplyr::select(-c(metric)) %>%
                          dplyr::filter(!is.na(group)) %>%
                          dplyr::select_if(~!(all(is.na(.)) | all(. == ""))) %>%
                          mutate(across(c(gene_symbol, group), factor)) %>%
                          distinct(.))) %>%
    bind_rows(., .id = "contrast")
  
  namevar <- deparse(substitute(rlist))
  assign(paste0(namevar, "_renamed_tpm_df"), df, envir = .GlobalEnv)
  rm(res_list)
  rm(df)
}


# Import data for DESeq analysis ------------------------------------------

##Prepare data
#import the preprocessed AnnotationHub reference

gene_tx_EnsDb <- read_rds(file = paste0(dir, "/data/gene_tx_EnsDb.rds"))

#import the counts files
counts_combined <- read_tsv(file = paste0(dir, "/data/counts_combined.tsv.gz"),
                            col_names = TRUE,
                            name_repair = "unique",
                            trim_ws = TRUE)

#import metadata
rnaseq_metadata <- read_tsv(file = paste0(dir, "/data/rnaseq_metadata.tsv"),
                            col_names = TRUE,
                            name_repair = "unique",
                            trim_ws = TRUE) %>%
  dplyr::mutate(treatment = factor(treatment, levels = c("DMSO", "Kw", "Taxol", "KwTaxol")),
                timepoint = factor(timepoint, levels = sort(unique(.$timepoint), decreasing = FALSE)),
                cell_type = factor(cell_type, levels = c("nociceptor", "brca")),
                cell_line = factor(cell_line, levels = c("iPSCs", "Hs578T", "MCF7", "MDAMB231"))) %>%
  dplyr::arrange(., by = treatment, timepoint, cell_type) %>%
  dplyr::mutate(sample_name = factor(sample_name, levels = unique(.$sample_name))) %>%
  dplyr::mutate(celltype_treatment_timepoint = factor(celltype_treatment_timepoint),
                celltype_treatment = factor(celltype_treatment),
                cellline_treatment_timepoint = factor(cellline_treatment_timepoint),
                SeqRun = factor(SeqRun))

#summarize transcripts to genes, remove non-protein-coding genes
#(multicore enabled)
counts_rep_protcod_df <- counts_combined %>%
  dplyr::mutate(gene_symbol = na_if(str_squish(gene_symbol), ""), 
                description = na_if(str_squish(description), "")) %>%
  dplyr::mutate(gene_symbol = factor(gene_symbol),
                gene_biotype = factor(gene_biotype, levels = levels(addNA(gene_tx_EnsDb$gene_biotype)), 
                                      labels = c(levels(gene_tx_EnsDb$gene_biotype), "Unspecified"), exclude = NULL)) %>%
  dplyr::slice(which(gene_symbol != "")) %>%
  dplyr::slice(which(gene_biotype == "protein_coding")) %>%
  droplevels(.) %>%
  tidyr::pivot_longer(., cols = !c(tx_id, gene_id, gene_symbol, description, gene_biotype),
                      names_to = "seq_name",
                      values_to = "reads") %>%
  dplyr::mutate(across(!contains("reads"), factor)) %>%
  dplyr::left_join(., rnaseq_metadata %>%
                     dplyr::select(seq_name, sample_name) %>%
                     distinct(.) %>%
                     droplevels(.) %>%
                     dplyr::mutate(across(everything(), factor))) %>%
  dplyr::select(-c(seq_name, tx_id)) %>%
  dplyr::mutate(across(!contains("reads"), factor)) %>%
  distinct(.) %>%
  multidplyr::partition(cluster) %>% #speed it up; but if you don't have multicore processing, comment this out
  dplyr::group_by(across(where(is.factor))) %>%
  dplyr::summarise(counts = sum(reads, na.rm = TRUE)) %>%
  collect() %>% #speed it up; see above
  droplevels(.) %>%
  as_tibble() %>%
  dplyr::mutate(sample_name = factor(sample_name, 
                                     levels = levels(rnaseq_metadata$sample_name)))

tpm_rep_protcod_df <- bind_rows(F_calculateTPM(cell = "brca", counts_df = counts_rep_protcod_df, mdata = rnaseq_metadata),
                                F_calculateTPM(cell = "nociceptor", counts_df = counts_rep_protcod_df, mdata = rnaseq_metadata))


# Run DESeq ---------------------------------------------------------------

##Run DESeq

deseq_design <- "~celltype_treatment_timepoint + SeqRun"
deseq_nociceptor <- F_deseq(cell = "nociceptor", 
                            counts_df = counts_rep_protcod_df, 
                            mdata = rnaseq_metadata,
                            deseq_design)

deseq_design <- "~cellline_treatment_timepoint + SeqRun"
deseq_brca <- F_deseq(cell = "brca", 
                      counts_df = counts_rep_protcod_df, 
                      mdata = rnaseq_metadata,
                      deseq_design)


# Extract DESeq results ---------------------------------------------------

##Extract results from DESeq

res_brca <- F_extract_results_multi(dlist = deseq_brca, 0.05)
res_nociceptor <- F_extract_results_multi(deseq_nociceptor, 0.05)


##Keep the results with p < 0.05
##these headers are derived from DESeq's default naming
headers <- c("contrast", "gene_symbol", "baseMean", "padj", "log2FoldChange_MMSE")
F_keep_sde(res_nociceptor, headers, 0.05)
F_keep_sde(res_brca, headers, 0.05)
#makes lists: 
#res_nociceptor_sig_list
#res_brca_sig_list


##Make your DESeq results legible.

F_rename_sde(deseq_nociceptor, res_nociceptor_sig_list)
F_rename_sde(deseq_brca, res_brca_sig_list)
#makes lists:
#res_nociceptor_sig_list_renamed_tpm
#res_brca_sig_list_renamed_tpm



# Save output from DESeq --------------------------------------------------

##Save all the files please
write_rds(res_nociceptor, paste0(dir, "/data/res_nociceptor.rds"), compress = "gz")
write_rds(res_brca, paste0(dir, "/data/res_brca.rds"), compress = "gz")
write_rds(res_nociceptor_sig_list, paste0(dir, "/data/res_nociceptor_sig_list.rds"), compress = "gz")
write_rds(res_brca_sig_list, paste0(dir, "/data/res_brca_sig_list.rds"), compress = "gz")
write_rds(res_nociceptor_sig_list_renamed_tpm, paste0(dir, "/data/res_nociceptor_sig_list_renamed_tpm.rds"), compress = "gz")
write_rds(res_brca_sig_list_renamed_tpm, paste0(dir, "/data/res_brca_sig_list_renamed_tpm.rds"), compress = "gz")


save(deseq_brca, res_brca, res_brca_sig_list, res_brca_sig_list_renamed_tpm, file = paste0(dir, "/data/res_brca_deseq_extra.RData"))
save(deseq_nociceptor, res_nociceptor, res_nociceptor_sig_list, res_nociceptor_sig_list_renamed_tpm, file = paste0(dir, "/data/res_nociceptor_deseq_extra.RData"))



