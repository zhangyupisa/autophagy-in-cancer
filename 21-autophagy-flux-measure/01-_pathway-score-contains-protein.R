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
pi3kakt_name <- c("AKT", "GSK", "p27", "PRAS40", "Tuberin_pT", "INPP4B", "PTEN")
pi3kakt_type <- c("p", "p", "p", "p", "p", "n", "n")
pi3kakt <- rep("PI3KAKT", length(pi3kakt_pat))

# RASMAPK
rasmapk_pat <- c("A.Raf_|c\\.Jun_pS73|C.Raf_|JNK_|MAPK_|MEK1_|p38_p|p90RSK_p|Shc_p|YB-1_p")
rasmapk_name <- c("cj")
rasmapk_type <- c("p")
rasmapk <- c("RASMAPK")

# RTK
rtk_pat <- c("HER2_|HER3_|Ret_|Shc_p|Src_p", "EGFR_")
rtk_name <- c("cj", "EGFR_")
rtk_type <- c("p", "p")
rtk <- c("RTK", "RTK")

# TSCmTOR
tscmtor_pat <- c("4E.BP1_p|mTOR_|p70S6K_|Rictor_", "S6_")
tscmtor_name <- c("cj", "pS6")
tscmtor_type <- c("p", "p")
tscmtor <- c("TSCmTOR", "TSCmTOR")

#Apoptosis 
apoptosis_pat <- c("Bak|Bid|Bim|Caspase.|Bax", "Bcl.2|Bcl.xL|Bad_|cIAP")
apoptosis_name <- c("cj", "n")
apoptosis_type <- c("p", "n")
apoptosis <- c("Apoptosis", "Apoptosis")

# Cell cycle
cellcycle_pat <- c("CDK1|Cyclin|p27_|PCNA")
cellcycle_name <- c("cj")
cellcycle_type <- c("p")
cellcycle <- c("CellCycle")

#DNA damage
dnad_pat <- c("53BP1|ATM|BRCA2|Chk1_|Chk2_|Ku80|Mre11|PARP|Rad50|Rad51|XRCC1|p53")
dnad_name <- c("cj")
dnad_type <- c("p")
dnad <- c("DNADamage")

# EMT 
emt_pat <- c("Collagen|Fibronectin|N.Cadherin", "Claudin.7|E.Cadherin")
emt_name <- c("cj", "cjl")
emt_type <- c("p", "n")
emt <- c("EMT", "EMT")

# ER
er_pat <- c("ER-|PR")
er_name <- c("cj")
er_type <- c("p")
er <- c("Hormone ER")

# AR
ar_pat <- c("AR|INPP4B|GATA3|Bcl.2")
ar_name <- c("cj")
ar_type <- c("p")
ar <- c("Hormone AR")


path_way <-
  tibble::tibble(
    pat = c(pi3kakt_pat, rasmapk_pat, rtk_pat, tscmtor_pat, apoptosis_pat, cellcycle_pat, dnad_pat, emt_pat, er_pat, ar_pat), 
    name = c(pi3kakt_name, rasmapk_name, rtk_name, tscmtor_name, apoptosis_name, cellcycle_name, dnad_name, emt_name, er_name, ar_name),
    type = c(pi3kakt_type, rasmapk_type, rtk_type, tscmtor_type, apoptosis_type, cellcycle_type, dnad_type, emt_type, er_type, ar_type),
    pathway = c(pi3kakt, rasmapk, rtk, tscmtor, apoptosis, cellcycle, dnad, emt, er, ar)
  )

tcpa_name %>% 
  dplyr::filter(stringr::str_detect(protein, pi3kakt_pat))
