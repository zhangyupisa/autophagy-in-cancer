
# Library -----------------------------------------------------------------

library(magrittr)
library(ggplot2)

# Path --------------------------------------------------------------------

path_project <- "/home/liucj/data/project/06-autophagy"
path_data <- file.path(path_project, "TCGA")
path_out <- file.path(path_project, "21-autophagy-flux")

# Load data ---------------------------------------------------------------

rppa_expr <- readr::read_rds(path = file.path(path_data, "pancan33-rppa-expr-v4-l4.rds.gz")) %>% 
  dplyr::filter(cancer_types != "PANCAN19") 
rppa_name <- readr::read_rds(path = file.path(path_data, "rppa-name-symbol.rds.gz"))
gene_list <- readr::read_rds(path = file.path(path_data, "rds-03-a-atg-lys-gene-list.rds.gz"))
pcc <- readr::read_tsv(file = file.path(path_data, "PanCanAtlasTumors_color_coded_by_organ_system_20170302.tsv"))


# Regulators --------------------------------------------------------------

# AKT pS473, AKT pT308, GSK3 pS9,  GSK3 pS21-pS9, PRAS40 pT246 TSC2 pT1462
PI3K_AKT <- "pS473|pT308|pS9|pS21S9|pT246|pT1462" %>% stringr::str_replace_all(pattern = "\\|", replacement = "$\\|") 

# mTOR pS2448, RICTOR pT1135, 4EBP1 pS65, 4EBP1 pT37T46, 4EBP1 pT70, S6K pT389, S6 pS235S236, S6 pS240S244
mTOR <- "pS2448|pT1135|pS65|pT37T46|pT70|pT389|pS235S236|pS240S244" %>% stringr::str_replace_all(pattern = "\\|", replacement = "$\\|")

rppa_expr %>% 
  dplyr::mutate(
    mtor = purrr::map(
      .x = expr,
      # mTor score
      .f = function(.x) {
        .x %>% 
          dplyr::filter(stringr::str_detect(protein, pattern = mTOR)) %>% 
          dplyr::select(-symbol) %>% 
          dplyr::mutate(
            direc = as.numeric(
              plyr::revalue(protein, c("MTOR_pS2448" = 1, "P70S6K_pT389" = 1, "RICTOR_pT1135" = 1, "S6_pS235S236" = 1, "S6_pS240S244" = 1, "X4EBP1_pS65" = 1, "X4EBP1_pT37T46" = 1, "X4EBP1_pT70" = 1)
            )
          )
        ) %>% 
          tidyr::gather(key = barcode, value = rppa, -protein, -direc) %>% 
          dplyr::group_by(barcode) %>% 
          dplyr::summarise(mtor_score = sum(rppa * direc))
      }
    )
  ) %>% 
  dplyr::mutate(
    pi3k_akt = purrr::map(
      .x = expr,
      # pi3k_akt score
      .f = function(.x) {
        .x %>% 
          dplyr::filter(stringr::str_detect(protein, pattern = PI3K_AKT)) %>% 
          dplyr::select(-symbol) %>% 
          dplyr::mutate(
            direc = as.numeric(
              plyr::revalue(protein, c("AKT_pS473" = 1, "AKT_pT308" = 1, "GSK3_pS9" = 1, "GSK3ALPHABETA_pS21S9" = 1, "PRAS40_pT246" = 1, "TUBERIN_pT1462" = 1))
            )
          ) %>% 
          tidyr::gather(key = barcode, value = rppa, -protein, -direc) %>% 
          dplyr::group_by(barcode) %>% 
          dplyr::summarise(pik_score = sum(rppa * direc))
      }
    )
  ) %>% 
  dplyr::select(-expr) %>% 
  dplyr::mutate(
    merge = purrr::map2(
      .x = mtor, 
      .y = pi3k_akt,
      .f = function(.x, .y) {
        .x %>% dplyr::inner_join(.y, by = "barcode")
      }
    )
  ) %>% 
  dplyr::select(1,4) %>% 
  tidyr::unnest() ->
  mtor_pi3k_score

# from cancer cell protein paper
glue::glue("{PI3K_AKT}$|{mTOR}$|P62|BECLIN|FOXO3A|RAPTOR|AMPKALPHA_pT172|AMPKALPHA|MTOR") -> PI3K_AKT_mTOR

rppa_expr %>% 
  dplyr::mutate(p62 = purrr::map(
    .x = expr,
    .f = function(.x) {
      .x %>% 
        dplyr::filter(
          stringr::str_detect(protein, pattern = PI3K_AKT_mTOR)
        ) %>% 
        dplyr::select(-symbol) %>% 
        tidyr::gather(key = barcode, value = rppa, -protein) %>% 
        tidyr::spread(key = protein, value = rppa)
    }
  )) %>% 
  dplyr::select(-expr) %>% 
  tidyr::unnest() %>% 
  dplyr::inner_join(mtor_pi3k_score, by = c("cancer_types", "barcode")) -> 
  p62_mtor_pi3k_score


# protein abundances ------------------------------------------------------

names(p62_mtor_pi3k_score)[-c(1,2)] -> proteins

fn_protein_abundences <- function(.n, .s) {
  
  p62_mtor_pi3k_score %>% 
    dplyr::select(1,2, .s) %>% 
    dplyr::rename(p = rlang::UQ(.s)) %>% 
    dplyr::mutate(
      p = dplyr::case_when(
        p > 3 ~ 3,
        p < -2 ~ -2,
        TRUE ~ p
      )
    ) -> .d
  
  lev <- c("TGCT", "LGG", "GBM", "PRAD", "UCS", "PAAD", "PCPG", "THCA", "KICH",
           "MESO", "CHOL", "THYM", "SARC", "KIRP", "COAD", "BLCA", "READ", "UCEC",
           "STAD", "CESC", "BRCA", "KIRC", "HNSC", "ESCA", "SKCM", "LUAD", "OV", 
           "LUSC", "UVM", "ACC", "LIHC", "DLBC") %>% 
    rev()
  
  # .d %>% 
  #   dplyr::group_by(cancer_types) %>% 
  #   dplyr::summarise(m = median(p)) %>% 
  #   dplyr::arrange(dplyr::desc(m)) %>% 
  #   dplyr::pull(cancer_types) -> lev
    
  .d %>% 
    dplyr::mutate(cancer_types = factor(cancer_types, lev)) %>% 
    ggplot(aes(x = cancer_types, y = p)) +
    stat_boxplot(geom = 'errorbar', width = 0.3) +
    geom_boxplot(outlier.colour = NA) +
    geom_point(
      aes(color = cancer_types), 
      position = position_jitter(width = 0.05), 
      alpha = 0.4,
      size = 0.8) +
    scale_color_manual(
      name = "Cancer Types",
      values = dplyr::slice(pcc, match(lev, cancer_types)) %>% dplyr::pull(color)
    ) +
    theme(
      axis.line = element_line(color = "black"),
      axis.text.x = element_text(hjust = 1, vjust = 0.5, angle = 90),
      
      panel.grid = element_blank(),
      panel.background = element_rect(fill = "white", color = NA),
      panel.spacing.x = unit(0, "lines")
    ) +
    guides(color = F) +
    labs(
      x = "", 
      y = glue::glue("RPPA protein aundances"), 
      title = glue::glue("{.s} protein expression")
    ) -> .p
  
  .dir <- file.path(path_out, "01-rppa-abundances")
  if (!dir.exists(.dir)) dir.create(.dir)
  
  ggsave(
    filename = glue::glue("fig-{.n}-{.s}-rppa-abundances.pdf"),
    plot = .p,
    device = "pdf",
    path = .dir,
    width = 7,
    height = 3
  )
}

proteins %>% 
  dplyr::as_tibble() %>% 
  tibble::rownames_to_column() %>% 
  dplyr::rename(.n = rowname, .s = value) %>% 
  purrr::pwalk(.f = fn_protein_abundences)


# p62 correlation ---------------------------------------------------------

fn_p62_corr <- function(.s) {
  p62_mtor_pi3k_score %>% 
    dplyr::group_by(cancer_types) %>% 
    dplyr::select(1, 2, rlang::UQ(.s), p62 = P62LCKLIGAND) %>% 
    dplyr::do(
      cor.test(
        dplyr::pull(., p62),
        dplyr::pull(., rlang::UQ(.s)),
        method = "spearman"
      ) %>% 
        broom::tidy()
    ) %>% 
    dplyr::ungroup() %>% 
    dplyr::select(1,coef = estimate, pval = p.value) %>% 
    dplyr::mutate(fdr = p.adjust(pval, method = "fdr")) %>% 
    tibble::add_column(vs = glue::glue("p62_vs_{.s}"), .before = 1)
}

proteins[-12] %>% 
  purrr::map(.f = fn_p62_corr) %>% 
  dplyr::bind_rows() %>% 
  dplyr::mutate(
    vs = stringr::str_replace(vs, "p62_vs_", "")
  ) %>% 
  ggplot(aes(x = coef, -log10(pval))) + 
  geom_point()
