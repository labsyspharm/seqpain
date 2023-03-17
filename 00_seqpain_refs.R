## Sharon Grim
## 2023

## 00_seqpain_refs.R

##This script prepares the reference databases for differential gene expression
##and gene set enrichment analyses.
##All scripts assume we are working in a local directory called "seqpain"
##and that data files are within "/seqpain/data/"
##


library(tidyverse)
library(BiocParallel)
library(BiocManager)
library(AnnotationDbi)
library(GenomicFeatures)
library(AnnotationHub)
library(ensembldb)


# Ensembl database from Annotation Hub ------------------------------------

#the RNAseq data were processed through Salmon against the EnsDb v.107

ah <- AnnotationHub()
## query(ah, c("EnsDb", "Homo sapiens", "107")) #AH104864
## drop <- c("miRNA", "lncRNA", "non_coding", "antisense", "rRNA", "rRNA_pseudogene", "macro_lncRNA", "3prime_overlapping_ncRNA", "bidirectional_promoter_lncRNA", "snoRNA", "scaRNA", "snRNA")
gene_tx_EnsDb <- full_join(
  (transcripts(ah[["AH104864"]], return.type = "data.frame") %>%
     dplyr::select(tx_id, gene_id, tx_id_version, tx_biotype)),
  (genes(ah[["AH104864"]], return.type = "data.frame") %>% 
     dplyr::select(gene_id, gene_name, description, gene_biotype))) %>%
  dplyr::mutate(gene_name = factor(gene_name),
         gene_biotype = factor(gene_biotype),
         tx_biotype = factor(tx_biotype))

rm(ah)

# Save output of Ensembl db processing ------------------------------------

write_rds(gene_tx_EnsDb, paste0(dir, "/data/gene_tx_EnsDb.rds"), compress = "gz")



# C2_CP database ----------------------------------------------------------

##2022/12 used this version: https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2022.1.Hs/c2.cp.v2022.1.Hs.json

MSigDB_C2_CP_pathways_v2022_1 <- rjson::fromJSON(file = "https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2022.1.Hs/c2.cp.v2022.1.Hs.json")

MSigDB_C2_CP_pathways_flat <- MSigDB_C2_CP_pathways_v2022_1 %>%
  sapply(., get, x = "geneSymbols")

MSigDB_C2_CP_pathways_renamer <- bind_cols(
  (MSigDB_C2_CP_pathways_v2022_1 %>% names),
  (MSigDB_C2_CP_pathways_v2022_1 %>%
     sapply(., get, x = "systematicName"))) %>%
  setNames(., c("pathway", "systematicName"))


MSigDB_C2_CP_pathways_v2022_genes <- MSigDB_C2_CP_pathways_v2022_1 %>%
  lapply(., get, x = "geneSymbols") %>%
  setNames(., MSigDB_C2_CP_pathways_v2022_1 %>% names) %>%
  map(., ~.x %>%
        tibble::enframe(., name = NULL, value = "geneSymbols")) %>%
  map2(., names(.),
       ~.x %>%
         dplyr::mutate(pathway = .y)) %>%
  imap_dfr(., ~.x %>%
             pivot_longer(cols = !c(pathway),
                          names_to = "metric",
                          values_to = "gene_symbol")) %>%
  dplyr::select(pathway, gene_symbol) %>%
  left_join(., MSigDB_C2_CP_pathways_renamer) %>%
  relocate(systematicName) %>%
  dplyr::mutate(across(everything(), factor)) %>%
  collect %>%
  droplevels


MSigDB_C2_CP_pathways_wide <- MSigDB_C2_CP_pathways_v2022_1 %>%
  lapply(., get, x = "geneSymbols") %>%
  setNames(., MSigDB_C2_CP_pathways_v2022_1 %>% names) %>%
  map(., ~.x %>%
        tibble::enframe(., name = NULL, value = "geneSymbols")) %>%
  map2(., names(.),
       ~.x %>%
         dplyr::mutate(pathway = .y)) %>%
  imap_dfr(., ~.x %>%
             pivot_longer(cols = !c(pathway),
                          names_to = "metric",
                          values_to = "gene_symbol")) %>%
  dplyr::select(pathway, gene_symbol) %>%
  left_join(., MSigDB_C2_CP_pathways_renamer) %>%
  dplyr::mutate(present = 1) %>%
  pivot_wider(., id_cols = c(pathway, systematicName), 
              names_from = "gene_symbol",
              values_fill = 0,
              values_from = "present") %>%
  column_to_rownames(., var = "systematicName") %>%
  dplyr::select(-pathway)


# Save output -------------------------------------------------------------


write_rds(MSigDB_C2_CP_pathways_flat, paste0(dir, "/data/MSigDB_C2_CP_pathways_flat.rds"), compress = "gz")
write.table(MSigDB_C2_CP_pathways_renamer, paste0(dir, "/data/MSigDB_C2_CP_pathways_renamer.tsv"), 
            sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

save(MSigDB_C2_CP_pathways_v2022_1, MSigDB_C2_CP_pathways_v2022_genes, MSigDB_C2_CP_pathways_renamer, MSigDB_C2_CP_pathways_wide, file = paste0(dir, "/data/MSigDB_C2_CP_extra.RData"))


