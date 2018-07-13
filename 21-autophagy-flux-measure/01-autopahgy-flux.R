
# Library -----------------------------------------------------------------

library(magrittr)
library(ggplot2)

# Path --------------------------------------------------------------------

path_project <- "/home/liucj/data/project/06-autophagy"
path_data <- file.path(path_project, "TCGA")
path_out <- file.path(path_project, "21-autophagy-flux")

# Load data ---------------------------------------------------------------
# load protein expression data
rppa_expr <- readr::read_rds(path = file.path(path_data, "pancan33-rppa-expr-v4-l4.rds.gz")) %>% 
  dplyr::filter(cancer_types != "PANCAN19") 

# load mrna expression pancan33-expr.rds.gz
mrna_expr <- readr::read_rds(path = file.path(path_data, "pancan33-expr.rds.gz"))

# protein name for the symbol
rppa_name <- readr::read_rds(path = file.path(path_data, "rppa-name-symbol.rds.gz"))

# autophagy gene list
gene_list <- readr::read_rds(path = file.path(path_data, "rds-03-a-atg-lys-gene-list.rds.gz"))

# pancan atlas cancer tissue color
pcc <- readr::read_tsv(file = file.path(path_data, "PanCanAtlasTumors_color_coded_by_organ_system_20170302.tsv"))

# rppa pathway score
pathway_score <- readr::read_rds(path = file.path(path_data, "pancan32-rppa-score.rds.gz"))

# rppa pathway gene list
pathway_gene_list <- readr::read_rds(path = file.path(path_data, "tcga-pathway-gene-list.rds.gz"))


# p62 correlates with SQSTM1 mRNA expression ------------------------------
# p62 rppa expression
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

# sqstm1 mrna expression
mrna_expr %>% 
  dplyr::mutate(
    sqstm1 = purrr::map(
      .x = expr,
      .f = function(.x) {
        .x %>% dplyr::select(-2) %>% dplyr::filter(symbol == "SQSTM1") 
      }
    )
  ) %>% 
  dplyr::select(-2) ->
  sqstm1_expr

# merge p62 protein and mrna expression by barcode
p62_rppa_expr %>% 
  dplyr::inner_join(sqstm1_expr, by = "cancer_types") %>% 
  dplyr::mutate(
    p62_sqstm1 = purrr::map2(
      .x = p62, 
      .y = sqstm1,
      .f = function(.x, .y) {
        .x %>% 
          tidyr::gather(key = "barcode", value = "p62", -protein) %>% 
          dplyr::select(-1) %>% 
          dplyr::mutate(barcode = substr(barcode, start = 1, stop = 16)) -> .xx
        .y %>% 
          tidyr::gather(key = "barcode", value = "sqstm1", -symbol) %>% 
          dplyr::select(-1) %>% 
          dplyr::mutate(barcode = substr(barcode, start = 1, stop = 16)) -> .yy
        .xx %>% dplyr::inner_join(.yy, by = "barcode")
        
      }
    )
  ) %>% 
  dplyr::select(-2, -3) %>% 
  dplyr::mutate(
    coef = purrr::map(
      .x = p62_sqstm1,
      .f = function(.x) {
        cor.test(.x$p62, .x$sqstm1, method = "pearson") %>% 
          broom::tidy() %>% 
          dplyr::select(coef = estimate, pval = p.value)
      }
    )
  ) %>% 
  tidyr::unnest(coef) ->
  p62_sqstm1_coef

p62_sqstm1_coef %>% dplyr::arrange(coef) %>% print(n = Inf)

# p62 correlates with pathway score ---------------------------------------
# p62 rppa expression correlates with the rppa score.

p62_rppa_expr %>% 
  dplyr::inner_join(pathway_score, by = "cancer_types") %>% 
  dplyr::mutate(
    corr = purrr::map2(
      .x = p62,
      .y = rppa,
      .f = function(.x, .y) {
        .x %>% 
          tidyr::gather(key = barcode, value = p62, -protein) %>% 
          dplyr::select(-protein) %>% 
          dplyr::mutate(barcode = substr(barcode, 1, 12)) %>% 
          dplyr::group_by(barcode) %>% 
          dplyr::summarise(p62 = mean(p62)) %>% 
          dplyr::ungroup() ->
          .p62
        .y %>% 
          tidyr::spread(key = pathway, value = score) %>% 
          dplyr::inner_join(.p62, by = "barcode") -> .d
        
        names(.d)[-1] -> .corrname
        # combination with 2 elements
        
        combn(x = .corrname, m = 2, simplify = F) -> .corrcomb
        .corrcomb %>% purrr::map(.f = function(.cn){
          cor.test(x = .d[[.cn[1]]], y = .d[[.cn[2]]], method = "pearson") %>% 
            broom::tidy() %>% 
            dplyr::select(coef = estimate, pval = p.value) %>% 
            tibble::add_column(vs = paste0(.cn, collapse = "#"), .before = 1) 
        }) %>% 
          dplyr::bind_rows() %>% 
          dplyr::filter(grepl(pattern = "p62", x = vs))
      }
    )
  ) %>% 
  dplyr::select(cancer_types, corr) %>% 
  tidyr::unnest() -> p62_pathway




