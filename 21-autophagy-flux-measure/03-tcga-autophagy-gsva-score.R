
# Load library ------------------------------------------------------------

library(magrittr)
library(ggplot2)

# Path --------------------------------------------------------------------

path_project <- "/home/liucj/data/project/06-autophagy"
path_data <- file.path(path_project, "TCGA")
path_out <- file.path(path_project, "21-autophagy-flux")



# load data ---------------------------------------------------------------
tcga_gene_symbol <- readr::read_rds(file.path(path_data, 'tcga-gene-symbol.rds.gz'))




# pancan33-expr.rds.gz
mrna <- readr::read_rds(path = file.path(path_data, "pancan33-expr.rds.gz"))

gene_list <- readr::read_rds(path = file.path(path_data, "rds-03-a-atg-lys-gene-list.rds.gz")) 


gene_list %>% 
  dplyr::select(symbol) %>% 
  dplyr::left_join(tcga_gene_symbol, by = 'symbol') ->
  gene_list_entrez

gene_sets <- list(
  atg_lys = dplyr::filter(gene_list, type %in% c("atg", "lys")),
  atg = dplyr::filter(gene_list, type == "atg"),
  lys = dplyr::filter(gene_list, type == "lys"),
  atg_core = dplyr::filter(gene_list, pathway == "autophagesome formation-core"),
  atg_form = dplyr::filter(gene_list, pathway == "autophagesome formation"),
  lys_com = dplyr::filter(gene_list, pathway == "lysosome composition"),
  lys_deg = dplyr::filter(gene_list, pathway == "lysosome degradation")
) %>% 
  purrr::map("symbol") %>% 
  purrr::map(
    .f = function(.x){
      gene_list_entrez %>% 
        dplyr::filter(symbol %in% .x) %>% 
        dplyr::pull(entrez_id) ->
        .v
      names(.v) <- .x
      .v
    }
  )
  



# GSVA --------------------------------------------------------------------

library(GSVA)

fn_gsva <- function(.x, .y, gene_sets = gene_sets) {
  # .x <- .te$cancer_types
  # .y <- .te$expr[[1]]
  # print(.x)

  # extract tumor
  names(.y)[c(-1, -2)] -> .b
  .b[substr(.b, start = 14, stop = 15) == "01"] -> .b

  .y %>%
    dplyr::select(2, .b) %>%
    dplyr::mutate(entrez_id = as.integer(entrez_id)) %>% 
    tidyr::drop_na() %>%
    dplyr::filter_at(
      .vars = dplyr::vars(dplyr::starts_with("T")), 
      .vars_predicate = dplyr::any_vars(. != 0)
    ) %>% 
    dplyr::mutate_at(
      .vars = dplyr::vars(dplyr::starts_with("T")),
      .funs = dplyr::funs(log2(. + 1))
    ) ->
    .d

  .d_mat <- as.matrix(.d[, -1])
  rownames(.d_mat) <- .d$entrez_id

  # the signature is 'matrix, list'
  # showMethods(f = 'gsva')
  # getMethod(f = 'gsva', signature = list(expr = 'matrix', gset.idx.list = 'list'))
  # for the gsva 'kcdf' the GAUSSION for microarray,
  # 'Poisson' for the rnaseq count,
  
  .es_dif <- gsva(
    expr = .d_mat, 
    gset.idx.list = gene_sets,
    kcdf = "Gaussian",
    method = "gsva",
    mx.diff = TRUE, 
    verbose = FALSE, 
    parallel.sz = 1
  )

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
  dplyr::mutate(
    gsva = purrr::map2(
      .x = cancer_types, 
      .y = expr, 
      .f = fn_gsva, 
      gene_sets = gene_sets
      )
    ) %>%
  dplyr::select(-expr) %>%
  dplyr::collect() %>%
  dplyr::ungroup() %>%
  dplyr::select(-PARTITION_ID) ->
gene_set_gsva

parallel::stopCluster(cl)

if (!file.exists(file.path(path_data, "tcga-autophagy-gsva.rds.gz"))) {
  readr::write_rds(x = gene_set_gsva, path = file.path(path_data, "tcga-autophagy-gsva.rds.gz"), compress = "gz")
}





# End ---------------------------------------------------------------------