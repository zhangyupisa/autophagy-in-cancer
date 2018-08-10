
# Load library ------------------------------------------------------------

library(magrittr)
library(ggplot2)

# Path --------------------------------------------------------------------

path_project <- "/home/liucj/data/project/06-autophagy"
path_data <- file.path(path_project, "TCGA")
path_out <- file.path(path_project, "21-autophagy-flux")

# load data ---------------------------------------------------------------
# pancan33-expr.rds.gz
mrna <- readr::read_rds(path = file.path(path_data, "pancan33-expr.rds.gz"))

gene_list <- readr::read_rds(path = file.path(path_data, "rds-03-a-atg-lys-gene-list.rds.gz"))

gene_sets <- list(
  atg_lys = dplyr::filter(gene_list, type %in% c("atg", "lys")),
  atg = dplyr::filter(gene_list, type == "atg"),
  lys = dplyr::filter(gene_list, type == "lys"),
  atg_core = dplyr::filter(gene_list, pathway == "autophagesome formation-core"),
  atg_form = dplyr::filter(gene_list, pathway == "autophagesome formation"),
  lys_com = dplyr::filter(gene_list, pathway == "lysosome composition"),
  lys_deg = dplyr::filter(gene_list, pathway == "lysosome degradation")
) %>% purrr::map("symbol")
#
# 



# GSVA --------------------------------------------------------------------

library(GSVA)

fn_gsva <- function(.x, .y, gene_sets = gene_sets){
  # .x <- .te$cancer_types
  # .y <- .te$filter_expr[[1]]
  print(.x)
  
  names(.y)[c(-1, -2)] -> .b
  .b[substr(.b, start = 14, stop = 14) == 0] -> .b
  
  .y %>% 
    tidyr::drop_na() %>% 
    dplyr::select(1, .b) %>% 
    dplyr::filter_at(.vars = dplyr::vars(dplyr::starts_with("T")), .vars_predicate = dplyr::any_vars(. != 0)) -> 
    .d
  
  .d_mat <- as.matrix(.d[,-1])
  rownames(.d_mat) <- .d$symbol
  
  .es_dif <- gsva(.d_mat, gene_sets, method = "gsva", mx.diff = TRUE, verbose = FALSE, parallel.sz = 1)
  
  .es_dif %>% 
    as.data.frame() %>% 
    tibble::as_tibble() %>% 
    tibble::add_column(set = rownames(.es_dif), .before = 1)
}


# mrna %>% 
#   dplyr::mutate(gsva = purrr::map2(.x = cancer_types, .y = expr, .f = fn_gsva, gene_sets = gene_sets)) %>% 
#   dplyr::select(-expr) -> 
#   gene_set_gsva


cl <- multidplyr::create_cluster(nrow(mrna))
mrna %>% 
  multidplyr::partition(cluster = cl) %>%
  multidplyr::cluster_library("magrittr") %>%
  multidplyr::cluster_library("GSVA") %>%
  multidplyr::cluster_assign_value("gene_sets", gene_sets) %>%
  multidplyr::cluster_assign_value("fn_gsva", fn_gsva) %>%
  dplyr::mutate(gsva = purrr::map2(.x = cancer_types, .y = expr, .f = fn_gsva, gene_sets = gene_sets)) %>%
  dplyr::select(-expr) %>% 
  dplyr::collect() %>%
  dplyr::ungroup() %>%
  dplyr::select(-PARTITION_ID) ->
  gene_set_gsva
  
parallel::stopCluster(cl)

if (!file.exists(file.path(path_data, "tcga-autophagy-gsva.rds.gz"))) 
  readr::write_rds(x = gene_set_gsva, path = file.path(path_data, "tcga-autophagy-gsva.rds.gz"), compress = "gz")





# End ---------------------------------------------------------------------


