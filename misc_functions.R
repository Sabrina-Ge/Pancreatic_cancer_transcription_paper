library(dplyr)

convert_human_to_mouse_homologues <- function(genelist){
  out <- readRDS("references/mouse_to_human_homologues_pdatranscriptomics.RDS")
  rownames(out) <- NULL
  out <- out %>% column_to_rownames("HGNC.symbol")
  return(out[genelist, "MGI.symbol"])
}

basal1 <- read.table("references/DE/basal1vsall_scoreslimma_DE.tsv", header = T, sep = "\t")
basal2 <- read.table("references/DE/basal2vsall_scoreslimma_DE.tsv", header = T, sep = "\t")
classical1 <- read.table("references/DE/classical1vsall_scoreslimma_DE.tsv", header = T, sep = "\t")
classical2 <- read.table("references/DE/classical2vsall_scoreslimma_DE.tsv", header = T, sep = "\t")

# Exclude genes highly expressed in stromal populations (Human Protein Atlas)
basal1_exclusion <- c('C9', 'CEL')
basal2_exclusion <- c('ORM1', 'FGB', 'LBP', 'FGG', 'F2', 'ALB', 'SLCO1B3', 'FGA', 'IGFBP1', 'ORM2', 'HP', 'UGT2B4', 'C8B', 'CRP', 'ASGR2', 'SLCO1B1', 'ITIH1', 'HPX', 'ITIH4')
classical1_exclusion <- c('CYP2C9', 'SPINK1', 'TM4SF4', 'NPC1L1', 'PRSS2', 'PRSS3', 'PRSS1', 'INS', 'REG1A')
classical2_exclusion <- c('REG1A', 'TTR', 'GP2', 'REG1B', 'INS', 'AKR1C4', 'SPINK1', 'PRSS2', 'PDX1', 'ANG')
exclusion_genes <- c(basal1_exclusion, basal2_exclusion, classical1_exclusion, classical2_exclusion)

create_topn_genelists <- function(top_n = 100, genelist_names = c("Basal1", "Basal2", "Classical1", "Classical2")){
  genelist <- list(basal1, basal2, classical1, classical2) %>%
    lapply(function(x)
      x %>%
        filter(adj.P.Val < 0.01 & logFC > 1) %>%
        mutate(pi_value = logFC * -1 * log10(adj.P.Val)) %>%
        arrange(desc(pi_value)) %>%
        pull(gene_name) %>% setdiff(exclusion_genes) %>% head(top_n)
    )
  
  names(genelist) <- genelist_names
  return(genelist)
}

