# List 10 pathway protein


# Library -----------------------------------------------------------------

library(magrittr)

# Path --------------------------------------------------------------------

path_project <- "/home/liucj/data/project/06-autophagy"
path_data <- file.path(path_project, "TCGA")

# Load data ---------------------------------------------------------------

rppa_name_symbol <- readr::read_rds(file.path(path_data, "rppa-name-symbol.rds.gz"))


# tcpa name ---------------------------------------------------------------

rppa_json_url <- "http://tcpaportal.org/tcpa/_design/basic/_show/annotation-antibody_list/annotation-antibody"

jsonlite::fromJSON(txt = rppa_json_url) %>% 
  dplyr::as_tibble() %>% 
  dplyr::select(`0`:`5`) %>% 
  dplyr::rename(
    protein = `0`,
    symbol = `1`,
    status = `2`,
    origin = `3`,
    source = `4`,
    catalog = `5`
  ) %>% 
  dplyr::mutate_all(.funs = dplyr::funs(stringr::str_trim)) %>% 
  dplyr::mutate(
    origin = stringr::str_sub(origin, end = 1),
    status = dplyr::recode(
      status,
      "Validated as ELISA" = "V",
      "Use with Caution" = "C",
      .default = "E"
    ),
    protein = stringr::str_replace(protein, pattern = "^x", "")
  ) %>% 
  dplyr::mutate(
    protein = stringr::str_c(protein, origin, status, sep = "-")
  ) %>% 
  dplyr::select(protein, symbol) %>% 
  dplyr::mutate(name = stringr::str_replace(protein, pattern = "-\\w-\\w$", "")) %>% 
  dplyr::mutate(name = stringr::str_to_upper(name)) %>% 
  dplyr::mutate(name = stringr::str_replace_all(name, pattern = "-|_", "")) -> tcpa_name

tcpa_name %>% readr::write_rds(path = file.path(path_data, "tcpa-rppa-name.rds.gz"), compress = "gz")


#pathways -------------------------------------------------------------------------------------
# pi3kakt
pi3kakt_pat <- c("Akt_", "GSK3", "p27", "PRAS40_", "Tuberin", "INPP4", "PTEN") %>% paste0(collapse = "|")

# RASMAPK
rasmapk_pat <- c("Raf|pS73|JNK|MAPK|MEK1|p38|p90RSK|Shc|YB")

# RTK
rtk_pat <- c("HER2|HER3|Ret|Shc|Src|^EGFR")


# TSCmTOR
tscmtor_pat <- c("^4E-BP1|mTOR|p70-S6K|Rictor|^S6")

#Apoptosis 
apoptosis_pat <- c("Bak|Bid|Bim|Caspase-|Bax|Bcl2|Bcl-xL|Bad|cIAP")

# Cell cycle
cellcycle_pat <- c("CDK1|Cyclin|p27|PCNA")

#DNA damage
dnad_pat <- c("53BP1|ATM|BRCA2|Chk1|Chk2|Ku80|Mre11|PARP|Rad50|Rad51|XRCC1|^p53")


# EMT 
emt_pat <- c("Collagen|Fibronectin|Cadherin|Claudin|E.Cadherin")


# ER
er_pat <- c("ER-|PR-")


# AR
ar_pat <- c("^AR-|INPP4B|GATA3|Bcl2")


# protein list ------------------------------------------------------------


path_way <-
  tibble::tibble(
    pattern = c(pi3kakt_pat, rasmapk_pat, rtk_pat, tscmtor_pat, apoptosis_pat, cellcycle_pat, dnad_pat, emt_pat, er_pat, ar_pat),
    pathway = c("PI3K/AKT", "RAS/MAPK", "RTK", "TSC/mTOR", "Apoptosis", "Cell Cycle", "DNA Damage Response", "EMT", "Hormone ER", "Hormone AR")
  )

path_way %>% 
  dplyr::mutate(
    gene_list = purrr::map(
      .x = pattern,
      .f = function(.x) {
        tcpa_name %>% 
          dplyr::filter(stringr::str_detect(string = protein, pattern = stringr::regex(.x, ignore_case = TRUE)))
      }
    )
  ) %>% 
  dplyr::select(2,3,1) -> 
  pathway_gene_list
pathway_gene_list %>% 
  readr::write_rds(path = file.path(path_data, "tcga-pathway-gene-list.rds.gz"), compress = "gz")

