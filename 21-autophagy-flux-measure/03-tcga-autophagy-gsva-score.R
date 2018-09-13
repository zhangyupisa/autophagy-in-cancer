
# Load library ------------------------------------------------------------

library(magrittr)
library(ggplot2)

# Path --------------------------------------------------------------------

path_project <- "/home/liucj/data/project/06-autophagy"
path_data <- file.path(path_project, "TCGA")
path_out <- file.path(path_project, "21-autophagy-flux")



# load data ---------------------------------------------------------------
tcga_gene_symbol <- readr::read_rds(file.path(path_data, 'tcga-gene-symbol.rds.gz'))

rppa_expr <- readr::read_rds(path = file.path(path_data, "pancan33-rppa-expr-v4-l4.rds.gz")) %>% 
  dplyr::filter(cancer_types != "PANCAN19")



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
  



# gsva score for tcga data ------------------------------------------------


mrna %>% 
  dplyr::mutate(
    tumor = purrr::map2_int(
      .x = cancer_types, 
      .y = expr,
      .f = function(.x, .y) {
        # .x <- .te$cancer_types
        # .y <- .te$expr[[1]]
        
        names(.y)[c(-1, -2)] -> .b
        .b[as.numeric(substr(.b, start = 14, stop = 15)) < 10] -> .b
        length(.b)
      }
    )
  ) %>% 
  print(n = Inf)

library(GSVA)

fn_gsva <- function(.x, .y, gene_sets = gene_sets) {
  # .x <- .te$cancer_types
  # .y <- .te$expr[[1]]
  # print(.x)

  # extract tumor
  names(.y)[c(-1, -2)] -> .b
  .b[as.numeric(substr(.b, start = 14, stop = 15)) < 10] -> .b

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


# gsva score distribution -------------------------------------------------



# p62 correlates with gsva score ------------------------------------------
rppa_expr %>% 
  dplyr::mutate(
    p62 = purrr::map(
      .x = expr, 
      .f = function(.x) {
        .x %>% dplyr::filter(protein == "P62LCKLIGAND") %>% dplyr::select(-1)
      })
  ) %>% 
  dplyr::select(-2) ->
  p62_rppa_expr


gene_set_gsva %>% 
  dplyr::inner_join(p62_rppa_expr, by = 'cancer_types') %>% 
  dplyr::mutate(
    p62_auto = purrr::map2(
      .x = p62,
      .y = gsva,
      .f = function(.x, .y) {
        .x %>% 
          tidyr::gather(key = "barcode", value = "p62", -protein) %>% 
          dplyr::select(-protein) %>% 
          dplyr::mutate(barcode = substr(x = barcode, start = 1, stop = 12)) %>% 
          dplyr::group_by(barcode) %>% 
          dplyr::summarise(p62 = mean(p62)) ->
          .xx
        
        .y %>% 
          # dplyr::filter(set == "atg_lys") %>% 
          tidyr::gather(key = "barcode", value = "score", -set) %>% 
          tidyr::spread(key = set, value = score) %>% 
          dplyr::filter(as.numeric(substr(x = barcode, start = 14, stop = 15)) < 10) %>% 
          dplyr::mutate(barcode = substr(x = barcode, start = 1, stop = 12)) %>% 
          dplyr::group_by(barcode) %>% 
          dplyr::summarise_if(.predicate = is.double, .funs = dplyr::funs(mean)) ->
          .yy
        
        .xx %>% dplyr::inner_join(.yy, by = "barcode")
      }
    )
  ) %>% 
  dplyr::select(1, 4) %>% 
  dplyr::mutate(
    corr = purrr::map(
      .x = p62_auto, 
      .f = function(.x) {
        names(gene_sets) %>% 
          purrr::map(
            .f = function(.z) {
              cor.test(.x[['p62']], .x[[.z]], method = "spearman", exact = FALSE) %>% 
                broom::tidy() %>% 
                dplyr::select(coef = estimate, pval = p.value)
            }
          ) ->
          .d_corr
        names(.d_corr) <- names(gene_sets)
        .d_corr %>% 
          tibble::enframe() %>% 
          tidyr::unnest()
      }
    )
  ) %>% 
  dplyr::select(1,3)  ->
  corr_p62_auto

corr_p62_auto %>%
  tidyr::unnest() %>% 
  dplyr::arrange(coef) %>% 
  dplyr::filter(pval < 0.05) %>% 
  print(n = Inf)

# End ---------------------------------------------------------------------
