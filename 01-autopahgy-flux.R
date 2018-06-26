
# Library -----------------------------------------------------------------

library(magrittr)

# Path --------------------------------------------------------------------

path_project <- "/home/liucj/data/project/06-autophagy"
path_data <- file.path(path_project, "TCGA")
path_out <- file.path(path_project, "21-autophagy-flux")

# Load data ---------------------------------------------------------------

rppa_expr_l4 <- readr::read_rds(path = file.path(path_data, "pancan33-rppa-expr-v4-l4.rds.gz"))
rppa_expr_l3 <- readr::read_rds(path = file.path(path_data, "pancan32-rppa-expr-v4-l3.rds.gz"))

rppa_expr_l3 %>% 
  dplyr::rename(expr_l3 = expr) %>% 
  dplyr::inner_join(rppa_expr_l4, by = "cancer_types") %>% 
  dplyr::rename(expr_l4 = expr) ->
  l3_l4

l3_l4 %>%
  dplyr::mutate(
    diff = purrr::map2(
      .x = expr_l3,
      .y = expr_l4,
      .f = function(.x, .y) {
        .x %>% 
          dplyr::filter(symbol == "SQSTM1") %>% 
          tidyr::gather(key = "barcode", value = "l3", -c(symbol, protein)) %>% 
          dplyr::select(3,4) ->
          .xx
        
        .y %>% 
          dplyr::filter(symbol == "P62LCKLIGAND") %>% 
          tidyr::gather(key = "barcode", value = "l4", -c(symbol, protein)) %>% 
          dplyr::select(3,4) ->
          .yy
        
        .xx %>% 
          dplyr::inner_join(.yy, by = "barcode") ->
          .zz
      }
    )
  )

# use l4

# p62 ---------------------------------------------------------------------

rppa_expr_l4 %>% 
  dplyr::filter(cancer_types != "PANCAN19") %>% 
  dplyr::mutate(
    p62 = purrr::map(
      .x = expr,
      .f = function(.x) {
        .x %>% 
          # P62LCKLIGAND is SQSTM1
          dplyr::filter(symbol == "P62LCKLIGAND")
      }
    )
  )
