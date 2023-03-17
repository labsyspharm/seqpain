## Sharon Grim
## 2023

## 02_seqpain_fgsea.R

##This script takes the prepared data from DESeq analyses
##and conducts fGSEA against the DMSO controls.
##All scripts assume we are working in a local directory called "seqpain"
##and that data files are within "/seqpain/data/"

library(tidyverse)
library(BiocParallel)
library(fgsea)

##if you have a modern computer, you can use multicore processing.
cluster <- multidplyr::new_cluster(4)
multidplyr::cluster_library(cluster, "dplyr")
bpparam_multi <- MulticoreParam(timeout = 1000,
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

#This retrieves all DESeq contrasts that were against DMSO (control) (multicore enabled)
F_get_vs_DMSO <- function(rlist, res_list, contrast_type, set.alpha){
  
  namevar <- deparse(substitute(rlist)) %>%
    gsub("(deseq_)|(res_)", "", .)
  tpm_avg_by_contrast <- res_list[["tpm_avg_by_contrast"]]
  renamer <- res_list[["renamer"]]
  res_sig_renamed_avgtpm <- res_list[["res_sig_renamed_avgtpm"]]
  
  log2FC_tpm <- left_join(
    (tpm_avg_by_contrast %>%
       ungroup() %>%
       dplyr::select(gene_symbol, avgTPM, 
                       # sdTPM, 
                     ({{contrast_type}})) %>%
       dplyr::filter(., !grepl(paste(keep, collapse="|"), {{contrast_type}}, perl = TRUE)) %>%
       dplyr::select_if(~!(all(is.na(.)) | all(. == ""))) %>%
       collect() %>%
       droplevels %>%
       dplyr::mutate(across(!contains("TPM"), factor)) %>%
       dplyr::rename(first = {{contrast_type}}) %>%
       collect() %>%
       group_by(gene_symbol, first) %>%
       left_join(., (renamer %>%
                       dplyr::slice(which(!grepl("_vs_", cell_line))) %>% #keep the contrasts wrt cell line
                       dplyr::select(label, renamed, first, second) %>%
                       dplyr::filter(., grepl(paste(keep, collapse="|"), second, perl = TRUE)) %>%
                       distinct(.keep_all = TRUE) %>%
                       droplevels))), 
    (tpm_avg_by_contrast %>%
       ungroup() %>%
       dplyr::select(gene_symbol, avgTPM, 
                     ({{contrast_type}})) %>%
       dplyr::filter(., grepl(paste(keep, collapse="|"), {{contrast_type}}, perl = TRUE)) %>%
       dplyr::select_if(~!(all(is.na(.)) | all(. == ""))) %>%
       collect() %>%
       droplevels %>%
       dplyr::mutate(across(!contains("TPM"), factor)) %>%
       dplyr::rename_with(., ~ gsub("TPM", "TPM.control", .x, fixed = TRUE),  contains("TPM")) %>%
       dplyr::rename(second = {{contrast_type}}) %>%
       distinct(.keep_all = TRUE))) %>%
    tidyr::drop_na(label, renamed) %>%
    collect() %>%
    group_by(label, renamed, gene_symbol) %>% 
    multidplyr::partition(cluster) %>%
    dplyr::mutate(across(ends_with("TPM"), ~ log2( (. + 1) / ((get(stringr::str_replace(cur_column(), 
                                                                                        "TPM$", "TPM.control"))) + 1)), .names = "log2FC_{.col}")) %>%
    collect() %>% 
    as_tibble() %>%
    dplyr::select(label, renamed, first, second, gene_symbol, avgTPM, avgTPM.control, log2FC_avgTPM) %>%
    pivot_longer(., cols = contains("TPM"),
                 names_to = c("parameter"),
                 values_to = "value") %>%
    relocate(c(label, renamed), .before = "first") %>%
    dplyr::mutate(renamed = gsub(":iPSCs", "", renamed)) %>%
    dplyr::mutate(renamed = gsub("_vs_DMSO", "", renamed)) %>%
    dplyr::mutate(renamed = gsub("_vs_0h", "_vs_DMSO:0h", renamed)) %>%
    dplyr::mutate(across(-c(value), factor)) %>%
    collect()
  
  sde_deseq <- res_sig_renamed_avgtpm %>%
    dplyr::filter(., grepl(paste(keep, collapse="|"), second, perl = TRUE)) %>%
    ungroup() %>%
    dplyr::select(label, renamed, gene_symbol, first, padj) %>%
    collect() %>%
    droplevels() %>%
    tidyr::drop_na() %>% 
    dplyr::mutate(renamed = gsub(":iPSCs", "", renamed)) %>%
    dplyr::mutate(renamed = gsub("_vs_DMSO", "", renamed)) %>%
    dplyr::mutate(renamed = gsub("_vs_0h", "_vs_DMSO:0h", renamed)) %>%
    dplyr::mutate(padj.sig = factor(ifelse(`padj` <= set.alpha, "sig", "insig"),
                                    levels = c("sig", "insig"))) %>%
    dplyr::mutate(padj.sig = factor(padj.sig, 
                                    levels = levels(addNA(padj.sig)),
                                    labels = c("sig", "insig", "cut"),
                                    exclude = NULL)) %>%
    dplyr::mutate(across(-c(padj), factor)) %>%
    collect()
  
  log2FC_deseq <- rlist %>%
    dplyr::filter(., contrast %in% levels(log2FC_tpm$label)) %>%
    ungroup() %>%
    dplyr::select(contrast, gene_symbol, log2FoldChange, log2FoldChange_MMSE) %>%
    collect() %>%
    droplevels() %>%
    tidyr::drop_na(log2FoldChange) %>%
    left_join(., (renamer %>%
                    dplyr::slice(which(!grepl("_vs_", cell_line))) %>% #keep the contrasts wrt cell line
                    dplyr::select(label, renamed, first, second) %>%
                    dplyr::filter(., grepl(paste(keep, collapse="|"), second, perl = TRUE)) %>%
                    distinct(.keep_all = TRUE) %>%
                    droplevels), by = c("contrast" = "label")) %>%
    dplyr::mutate(renamed = gsub(":iPSCs", "", renamed)) %>%
    dplyr::mutate(renamed = gsub("_vs_DMSO", "", renamed)) %>%
    dplyr::mutate(renamed = gsub("_vs_0h", "_vs_DMSO:0h", renamed)) %>%
    collect() %>%
    pivot_longer(., cols = contains("log2FoldChange"),
                 names_to = c("parameter"),
                 values_to = "value") %>%
    dplyr::rename(label = contrast) %>%
    dplyr::mutate(across(-c(value), factor)) %>%
    collect()
  
  output_list <- list(sde_deseq, log2FC_tpm, log2FC_deseq) %>%
    setNames(., c("sde_deseq", "log2FC_tpm", "log2FC_deseq"))
  assign(paste0("res_", namevar, "_log2FC_treatment_list"), output_list, envir = .GlobalEnv)
  
  rm(sde_deseq)
  rm(log2FC_tpm)
  rm(log2FC_deseq)
}


#This performs fast GSEA (multicore enabled)
F_gsea_list <- function(resultslist, comparisons, renamelist, calc1, calc2){
  calc1 <- rlang::parse_expr(calc1)
  calc2 <- rlang::parse_expr(calc2)
  df_of_sde <- resultslist[["sde_deseq"]]
  df_of_calc1 <- resultslist[[{calc1}]]
  df_of_calc2 <- resultslist[[{calc2}]]
  
  temp_list2 <- df_of_sde %>%
    dplyr::select(gene_symbol, label) %>%
    droplevels
  temp_list1 <- df_of_calc1 %>%
    dplyr::filter(., grepl("log2", parameter)) %>%
    droplevels %>%
    dplyr::select(value, gene_symbol, label) %>%
    dplyr::filter(!is.na(value)) %>%
    dplyr::filter(label %in% {{comparisons}}) %>%
    collect %>%
    droplevels %>%
    dplyr::filter(!is.na(value)) %>%
    ungroup %>%
    arrange(by = dplyr::desc(value)) %>%
    as_tibble
  temp_list3 <- df_of_calc2 %>%
    dplyr::filter(., grepl("MMSE", parameter)) %>%
    droplevels %>%
    dplyr::select(value, gene_symbol, label) %>%
    dplyr::filter(!is.na(value)) %>%
    dplyr::filter(label %in% {{comparisons}}) %>%
    collect %>%
    droplevels %>%
    dplyr::filter(!is.na(value)) %>%
    ungroup %>%
    arrange(by = dplyr::desc(value)) %>%
    as_tibble
  temp_list4 <- temp_list2 %>%
    inner_join(., temp_list3) %>%
    dplyr::select(value, gene_symbol, label) %>%
    droplevels %>%
    as_tibble
  temp_list2 <- temp_list2 %>%
    inner_join(., temp_list1) %>%
    dplyr::select(value, gene_symbol, label) %>%
    droplevels %>%
    as_tibble
  
  renamer <- renamelist %>%
    dplyr::select(label, renamed) %>%
    distinct(label, .keep_all = TRUE) %>%
    dplyr::filter(., (.$label %in% levels(temp_list1$label))) %>%
    dplyr::select(label) %>%
    distinct %>%
    droplevels %>%
    unlist
  
  
  res1 <- map(renamer %>% set_names(),
              ~(dplyr::filter(temp_list1, grepl(.x, label)) %>%
                  dplyr::select(value, gene_symbol))) %>%
    map(., with, set_names(value, gene_symbol)) %>%
    map(., unlist, use.names = TRUE) %>%
    map(., ~fgsea::fgseaMultilevel(pathways=MSigDB_C2_CP_pathways_flat,
                                   stats=.x,
                                   nPermSimple=1e4,
                                   # nproc=1)
                                   BPPARAM=bpparam_multi) #multicore enabled
    ) %>%
    map(as_tibble)
  res2 <- map(renamer %>% set_names(),
              ~(dplyr::filter(temp_list2, grepl(.x, label)) %>%
                  dplyr::select(value, gene_symbol))) %>%
    map(., with, set_names(value, gene_symbol)) %>%
    map(., unlist, use.names = TRUE) %>%
    map(., ~fgsea::fgseaMultilevel(pathways=MSigDB_C2_CP_pathways_flat,
                                   stats=.x,
                                   nPermSimple=1e4,
                                   # nproc=1)
                                   BPPARAM=bpparam_multi) #multicore enabled
    ) %>%
    map(as_tibble)
  res3 <- map(renamer %>% set_names(),
              ~(dplyr::filter(temp_list3, grepl(.x, label)) %>%
                  dplyr::select(value, gene_symbol))) %>%
    map(., with, set_names(value, gene_symbol)) %>%
    map(., unlist, use.names = TRUE) %>%
    map(., ~fgsea::fgseaMultilevel(pathways=MSigDB_C2_CP_pathways_flat,
                                   stats=.x,
                                   nPermSimple=1e4,
                                   # nproc=1)
                                   BPPARAM=bpparam_multi) #multicore enabled
    ) %>%
    map(as_tibble)
  res4 <- map(renamer %>% set_names(),
              ~(dplyr::filter(temp_list4, grepl(.x, label)) %>%
                  dplyr::select(value, gene_symbol))) %>%
    map(., with, set_names(value, gene_symbol)) %>%
    map(., unlist, use.names = TRUE) %>%
    map(., ~fgsea::fgseaMultilevel(pathways=MSigDB_C2_CP_pathways_flat,
                                   stats=.x,
                                   nPermSimple=1e4,
                                   # nproc=1)
                                   BPPARAM=bpparam_multi) #multicore enabled
    ) %>%
    map(as_tibble)
  fgsea_namer <- c(paste0("all_", {calc1}),
                   paste0("sde_", {calc1}),
                   paste0("all_", {calc2}),
                   paste0("sde_", {calc2}))
  
  temp_list <- list(res1, res2, res3, res4) %>%
    setNames(., fgsea_namer)
  namevar <- deparse(substitute(resultslist)) %>%
    gsub("(^[[:alpha:]]+\\_[[:alpha:]]+).*", "\\1", ., perl = TRUE)
  
  assign(paste0(namevar, "_fgsea"), temp_list, envir = .GlobalEnv)
  
}

#This extracts only the significant results from GSEA
F_gsea_sig_paths <- function(fgsea, renamelist, set.alpha){
  map2(fgsea, fgsea %>% names,
       ~imap_dfr(., ~ .x %>%
                   dplyr::mutate(label = .y))) %>%
    map(., dplyr::select, label, pathway, pval, padj, NES, size) %>%
    imap_dfr(., ~ .x %>%
               mutate(subset = .y) %>%
               pivot_longer(cols = !c(label, subset, pathway),
                            names_to = "metric",
                            values_to = "value")) %>%
    dplyr::group_by(label, subset, pathway, metric) %>%
    pivot_wider(names_from = metric, values_from = value,
                id_cols = c(label, subset, pathway)) %>%
    dplyr::arrange(by = desc(-padj)) %>%
    dplyr::filter(padj <= set.alpha) %>%
    ungroup %>%
    dplyr::select(label, subset, pathway, padj, NES, size) %>%
    left_join(., renamelist %>%
                dplyr::select(label, renamed) %>%
                distinct(label, .keep_all = TRUE) %>%
                droplevels) %>%
    collect() %>%
    droplevels %>%
    distinct(.keep_all = TRUE)  %>%
    relocate(renamed) %>%
    dplyr::mutate(across(-c(padj, NES, size), factor)) %>%
    split(., f = .$subset, drop = TRUE) %>%
    map(., ~.x %>%
          collect)
}

#This gets a full matrix of NES and padj of pathways that were found at least once to be sig in a contrast
F_gsea_sig_wide <- function(fgsea, renamelist, set.alpha){
  
  #get a list of the pathways that were found sig within each subset
  temp_list <- map2(fgsea, fgsea %>% names,
                    ~imap_dfr(., ~ .x %>%
                                dplyr::mutate(label = .y))) %>%
    map(., dplyr::select, label, pathway, padj) %>%
    imap_dfr(., ~ .x %>%
               mutate(subset = .y) %>%
               pivot_longer(cols = !c(subset, label, pathway),
                            names_to = "metric",
                            values_to = "padj")) %>%
    dplyr::arrange(by = desc(-padj)) %>%
    dplyr::filter(padj <= set.alpha) %>%
    ungroup %>%
    collect %>%
    dplyr::select(subset, pathway) %>%
    distinct(.keep_all = TRUE) %>%
    dplyr::mutate(across(everything(), factor)) %>%
    split(., f = .$subset, drop = TRUE) %>%
    map(., ~.x %>%
          droplevels)
  
  #list of NES and padj of pathways in contrasts in each subset
  temp_list2 <- map2(fgsea, fgsea %>% names,
                     ~imap_dfr(., ~ .x %>%
                                 dplyr::mutate(label = .y))) %>%
    map(., dplyr::select, label, pathway, padj, NES) %>%
    imap_dfr(., ~ .x %>%
               mutate(subset = .y) %>%
               pivot_longer(cols = !c(label, subset, pathway),
                            names_to = "metric",
                            values_to = "value")) %>%
    dplyr::filter(subset %in% names(temp_list)) %>%
    collect %>%
    droplevels %>%
    dplyr::group_by(label, subset, pathway, metric) %>%
    pivot_wider(names_from = metric, values_from = value,
                id_cols = c(label, subset, pathway)) %>%
    ungroup %>%
    dplyr::select(label, subset, pathway, padj, NES) %>%
    left_join(., renamelist %>%
                dplyr::select(label, renamed) %>%
                distinct(label, .keep_all = TRUE) %>%
                droplevels) %>%
    collect() %>%
    droplevels %>%
    distinct(.keep_all = TRUE)  %>%
    relocate(renamed) %>%
    dplyr::mutate(across(-c(padj, NES), factor)) %>%
    split(., f = .$subset, drop = TRUE) %>%
    map(., ~.x %>%
          droplevels)
  
  temp_list3 <- map(temp_list %>% names,
                    ~temp_list2[[.x]] %>%
                      as_tibble %>%
                      dplyr::filter(pathway %in% temp_list[[.x]]$pathway) %>%
                      collect %>%
                      droplevels) %>%
    setNames(., temp_list %>% names)
  
  namevar <- deparse(substitute(fgsea))
  assign(paste0(namevar, "_sig_wide"), temp_list3, envir = .GlobalEnv)
  rm(temp_list)
  rm(temp_list2)
  rm(temp_list3)
  
}



# Import your reference databases -----------------------------------------

MSigDB_C2_CP_pathways_flat <- read_rds(file = paste0(dir, "/data/MSigDB_C2_CP_pathways_flat.rds"))


# Import your DESeq results -----------------------------------------------

res_nociceptor <- read_rds(file = paste0(dir, "/data/res_nociceptor.rds"))
res_nociceptor_sig_list_renamed_tpm <- read_rds(file = paste0(dir, "/data/res_nociceptor_sig_list_renamed_tpm.rds"))

res_brca <- read_rds(file = paste0(dir, "/data/res_brca.rds"))
res_brca_sig_list_renamed_tpm <- read_rds(file = paste0(dir, "/data/res_brca_sig_list_renamed_tpm.rds"))




# Prepare the data for fGSEA ----------------------------------------------

#prepare lists of gene expressions for comparisons against control (DMSO @ 0h or the corresponding timepoint)

keep <- c("2A", "2K", "2L", "2M")
F_get_vs_DMSO(res_nociceptor, res_nociceptor_sig_list_renamed_tpm, celltype_treatment_timepoint, 0.05)

write_rds(res_nociceptor_log2FC_treatment_list, paste0(dir, "/data/res_nociceptor_log2FC_treatment_list.rds"), compress = "gz")


F_get_vs_DMSO(res_brca, res_brca_sig_list_renamed_tpm, cellline_treatment_timepoint, 0.05)

write_rds(res_brca_log2FC_treatment_list, paste0(dir, "/data/res_brca_log2FC_treatment_list.rds"), compress = "gz")


# Conduct gene set enrichment analysis using fGSEA ------------------------
##This can be computationally intensive. 
##I recommend running an RStudio interactive session on the O2 compute cluster
##to perform these GSEA functions. I had success with 4 cores, 32Gb RAM,
##and 4 hours wall time
##Because parallelization is different on the compute cluster,
##the O2 versions of these functions are in a separate script.

keep <- c("2A", "2K", "2L", "2M")
idx <- res_nociceptor_sig_list_renamed_tpm[["renamer"]] %>% #looking at comparisons between timepoints and treatments within same cell line
  dplyr::filter(., grepl(paste(keep, collapse="|"), second, perl = TRUE)) %>%
  dplyr::select(label, renamed, cell_line, treatment, time) %>%
  dplyr::filter(., !grepl("_vs_", cell_line) & grepl("_vs_DMSO", treatment)) %>%
  droplevels %>%
  dplyr::mutate(across(.cols = everything(), as.factor)) %>%
  dplyr::select(label) %>%
  simplify %>%
  as.character

F_gsea_list(resultslist = res_nociceptor_log2FC_treatment_list,
            comparisons = idx,
            renamelist = res_nociceptor_sig_list_renamed_tpm[["renamer"]],
            calc1 = "log2FC_tpm",
            calc2 = "log2FC_deseq")
##res_nociceptor_fgsea

#having issues with brca because two samples are duplicates of each other :(
{
keep <- c("2A", "2K", "2L", "2M")
# idx <- res_brca_sig_list_renamed_tpm[["renamer"]] %>% #looking at comparisons between timepoints and treatments within same cell line
#   dplyr::filter(., grepl(paste(keep, collapse="|"), second, perl = TRUE)) %>%
#   dplyr::select(label, renamed, cell_line, treatment, time) %>%
#   dplyr::filter(., !grepl("_vs_", cell_line)) %>%
#   dplyr::filter(if_any(c(treatment, time), ~ !grepl("_vs_", .x))) %>%
#   droplevels %>%
#   dplyr::mutate(across(.cols = everything(), as.factor)) %>%
#   dplyr::select(label) %>%
#   simplify %>%
#   as.character

rm(res_brca_fgsea)
rm(temp_list)

keep <- c("2L")
idx <- res_brca_sig_list_renamed_tpm[["renamer"]] %>% #looking at comparisons between timepoints and treatments within same cell line
  dplyr::filter(., grepl(paste(keep, collapse="|"), second, perl = TRUE)) %>%
  dplyr::select(label, renamed, cell_line, treatment, time) %>%
  dplyr::filter(., !grepl("_vs_", cell_line)) %>%
  dplyr::filter(if_any(c(treatment, time), ~ !grepl("_vs_", .x))) %>%
  droplevels %>%
  dplyr::mutate(across(.cols = everything(), as.factor)) %>%
  dplyr::select(label) %>%
  simplify %>%
  as.character

temp_list <- F_gsea_list(resultslist = res_brca_log2FC_treatment_list,
            comparisons = idx,
            renamelist = res_brca_sig_list_renamed_tpm[["renamer"]],
            calc1 = "log2FC_tpm",
            calc2 = "log2FC_deseq")

keep <- c("2M")
idx <- res_brca_sig_list_renamed_tpm[["renamer"]] %>% #looking at comparisons between timepoints and treatments within same cell line
  dplyr::filter(., grepl(paste(keep, collapse="|"), second, perl = TRUE)) %>%
  dplyr::select(label, renamed, cell_line, treatment, time) %>%
  dplyr::filter(., !grepl("_vs_", cell_line)) %>%
  dplyr::filter(if_any(c(treatment, time), ~ !grepl("_vs_", .x))) %>%
  droplevels %>%
  dplyr::mutate(across(.cols = everything(), as.factor)) %>%
  dplyr::select(label) %>%
  simplify %>%
  as.character

temp_list2 <- F_gsea_list(resultslist = res_brca_log2FC_treatment_list,
                         comparisons = idx,
                         renamelist = res_brca_sig_list_renamed_tpm[["renamer"]],
                         calc1 = "log2FC_tpm",
                         calc2 = "log2FC_deseq")

}

{#2K seems to be the issue...
keep <- c("2K")
idx <- res_brca_sig_list_renamed_tpm[["renamer"]] %>% #looking at comparisons between timepoints and treatments within same cell line
  dplyr::filter(., grepl(paste(keep, collapse="|"), second, perl = TRUE)) %>%
  dplyr::select(label, renamed, cell_line, treatment, time) %>%
  dplyr::filter(., !grepl("_vs_", cell_line) & grepl("_vs_DMSO", treatment)) %>%
  droplevels %>%
  dplyr::mutate(across(.cols = everything(), as.factor)) %>%
  dplyr::select(label) %>%
  simplify %>%
  as.character

temp_list2 <- res_brca_log2FC_treatment_list[["sde_deseq"]] %>%
  dplyr::select(gene_symbol, label) %>%
  droplevels
temp_list1 <- res_brca_log2FC_treatment_list[["log2FC_tpm"]] %>%
  dplyr::filter(., grepl("log2", parameter)) %>%
  droplevels %>%
  dplyr::select(value, gene_symbol, label) %>%
  dplyr::filter(!is.na(value)) %>%
  dplyr::filter(label %in% idx) %>%
  collect %>%
  droplevels %>%
  dplyr::filter(!is.na(value)) %>%
  ungroup %>%
  arrange(by = dplyr::desc(value)) %>%
  as_tibble
temp_list3 <- res_brca_log2FC_treatment_list[["log2FC_deseq"]] %>%
  dplyr::filter(., grepl("MMSE", parameter)) %>%
  droplevels %>%
  dplyr::select(value, gene_symbol, label) %>%
  dplyr::filter(!is.na(value)) %>%
  dplyr::filter(label %in% idx) %>%
  collect %>%
  droplevels %>%
  dplyr::filter(!is.na(value)) %>%
  ungroup %>%
  arrange(by = dplyr::desc(value)) %>%
  as_tibble
temp_list4 <- temp_list2 %>%
  inner_join(., temp_list3) %>%
  dplyr::select(value, gene_symbol, label) %>%
  droplevels %>%
  as_tibble
temp_list2 <- temp_list2 %>%
  inner_join(., temp_list1) %>%
  dplyr::select(value, gene_symbol, label) %>%
  droplevels %>%
  as_tibble

#temp_list1: all 27 contrasts with 2K* as the second
#temp_list2: 10
#temp_list3: 15
#temp_list4: 10

g <- ggplot(data = temp_list1)
g <- g + theme_bw()
g <- g + geom_violin(aes(y = value, x = label),
                     draw_quantiles = c(0.5),
                     alpha = 0.3, position = position_dodge(0.8))
g
#cellline_treatment_timepoint_2S3_vs_2K3 seems to be the problem.
df <- temp_list1 %>%
  group_by(label) %>%
  dplyr::summarise(avg_value = mean(value, na.rm = TRUE))

df <- temp_list1 %>%
  dplyr::filter(label == "cellline_treatment_timepoint_2S3_vs_2K3") %>%
  droplevels




keep <- c("2K")
df <- res_brca_sig_list_renamed_tpm[["renamer"]] %>% #looking at comparisons between timepoints and treatments within same cell line
  dplyr::filter(., grepl(paste(keep, collapse="|"), second, perl = TRUE)) %>%
  dplyr::select(label, renamed, cell_line, treatment, time) %>%
  dplyr::filter(., !grepl("_vs_", cell_line) & grepl("_vs_DMSO", treatment)) %>%
  droplevels

renamer <- res_brca_sig_list_renamed_tpm[["renamer"]] %>%
  dplyr::select(label, renamed) %>%
  distinct(label, .keep_all = TRUE) %>%
  dplyr::filter(., (.$label %in% levels(temp_list1$label))) %>%
  dplyr::select(label) %>%
  distinct %>%
  droplevels %>%
  unlist


res1 <- map(renamer %>% set_names(),
            ~(dplyr::filter(temp_list1, grepl(.x, label)) %>%
                dplyr::select(value, gene_symbol))) %>%
  map(., with, set_names(value, gene_symbol)) %>%
  map(., unlist, use.names = TRUE) %>%
  map(., ~fgsea::fgseaMultilevel(pathways=MSigDB_C2_CP_pathways_flat,
                                 stats=.x,
                                 nPermSimple=1e4,
                                 # nproc=1)
                                 BPPARAM=bpparam_multi) #multicore enabled
  ) %>%
  map(as_tibble)


}

idx <- res_brca_sig_list_renamed_tpm[["renamer"]] %>% #looking at comparisons between timepoints and treatments within same cell line
  dplyr::filter(., grepl(paste(keep, collapse="|"), second, perl = TRUE)) %>%
  dplyr::select(label, renamed, cell_line, treatment, time) %>%
  dplyr::filter(., !grepl("_vs_", cell_line)) %>%
  dplyr::filter(if_any(c(treatment, time), ~ !grepl("_vs_", .x))) %>%
  droplevels %>%
  dplyr::mutate(across(.cols = everything(), as.factor)) %>%
  dplyr::select(label) %>%
  simplify %>%
  as.character

F_gsea_list(resultslist = res_brca_log2FC_treatment_list,
                          comparisons = idx,
                          renamelist = res_brca_sig_list_renamed_tpm[["renamer"]],
                          calc1 = "log2FC_tpm",
                          calc2 = "log2FC_deseq")

write_rds(res_brca_fgsea, paste0(dir, "/data/res_brca_fgsea.rds"), compress = "gz")

##res_brca_fgsea

# Retain only significant results and process for figures etc -------------

res_nociceptor_fgsea_sig <- F_gsea_sig_paths(res_nociceptor_fgsea, res_nociceptor_sig_list_renamed_tpm[["renamer"]], 0.05)

F_gsea_sig_wide(res_nociceptor_fgsea, res_nociceptor_sig_list_renamed_tpm[["renamer"]], 0.05)
##res_nociceptor_fgsea_sig_wide


res_brca_fgsea_sig <- F_gsea_sig_paths(res_brca_fgsea, res_brca_sig_list_renamed_tpm[["renamer"]], 0.05)

F_gsea_sig_wide(res_brca_fgsea, res_brca_sig_list_renamed_tpm[["renamer"]], 0.05)
##res_brca_fgsea_sig_wide



# Save output from fast GSEA ----------------------------------------------

write_rds(res_nociceptor_fgsea_sig_wide, paste0(dir, "/data/res_nociceptor_fgsea_sig_wide.rds"), compress = "gz")
write_rds(res_nociceptor_fgsea_sig, paste0(dir, "/data/res_nociceptor_fgsea_sig.rds"), compress = "gz")


write_rds(res_brca_fgsea_sig_wide, paste0(dir, "/data/res_brca_fgsea_sig_wide.rds"), compress = "gz")
write_rds(res_brca_fgsea_sig, paste0(dir, "/data/res_brca_fgsea_sig.rds"), compress = "gz")


save(res_brca_fgsea, res_brca_fgsea_sig_wide, res_brca_fgsea_sig, file = paste0(dir, "/data/res_brca_fgsea_extra.RData"))
save(res_nociceptor_fgsea, res_nociceptor_fgsea_sig_wide, res_nociceptor_fgsea_sig, file = paste0(dir, "/data/res_nociceptor_fgsea_extra.RData"))



