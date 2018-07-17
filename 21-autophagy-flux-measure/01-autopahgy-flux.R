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
pcc <- readr::read_tsv(file = file.path(path_data, "PanCanAtlasTumors-color-coded-by-organ-system-20170302.tsv"))

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
        cor.test(.x$p62, log2(.x$sqstm1 + 1), method = "spearman", exact = FALSE) %>% 
          broom::tidy() %>% 
          dplyr::select(coef = estimate, pval = p.value)
      }
    )
  ) %>% 
  tidyr::unnest(coef) ->
  p62_sqstm1_coef

human_read <- function(.x){
  if (.x > 0.1) {
    .x %>% signif(digits = 2) %>% toString()
  } else if (.x < 0.1 && .x > 0.001 ) {
    .x %>% signif(digits = 1) %>% toString()
  } else {
    .x %>% format(digits = 2, scientific = TRUE)
  }
}

p62_sqstm1_coef %>% 
  dplyr::select(1, 2) %>% 
  tidyr::unnest() %>% 
  dplyr::mutate(sqstm1 = log2(sqstm1 + 1)) -> 
  p62_sqstm1_expr

cor.test(p62_sqstm1_expr$p62, p62_sqstm1_expr$sqstm1, method = "spearman", exact = FALSE) %>% 
  broom::tidy() %>% 
  dplyr::select(ceof = estimate, pval = p.value) %>% 
  unlist() -> 
  p62_sqstm1_total_corr


p62_sqstm1_coef %>% 
  dplyr::select(-2) %>% 
  dplyr::inner_join(pcc, by = "cancer_types") %>% 
  dplyr::mutate(label = glue::glue("{cancer_types} R={format(coef, digits = 2)}")) %>% 
  dplyr::select(1, 2, 3, 4, 6) %>% 
  dplyr::mutate(pval = purrr::map_chr(pval, human_read)) %>% 
  dplyr::mutate(
    label = purrr::map2(
      .x = pval,
      .y = label,
      .f = function(.x, .y) {
        if (grepl(pattern = "e", x = .x)) {
          sub("-0", "-", strsplit(split = "e", x = .x, fixed = TRUE)[[1]]) -> .xx
          latex2exp::TeX(glue::glue("<<.y>> p=$<<.xx[1]>> \\times 10^{<<.xx[2]>>}$", .open = "<<", .close = ">>"))
        } else {
          latex2exp::TeX(glue::glue("{.y} p={.x}"))
        }
      }
    )
  ) %>% 
  dplyr::arrange(dplyr::desc(coef)) ->
  p62_sqstm1_coef_label


p62_sqstm1_expr %>% 
  dplyr::mutate(sqstm1 = ifelse(sqstm1 > 18.6, 18.6, sqstm1)) %>%
  dplyr::mutate(sqstm1 = ifelse(sqstm1 < 9, 9, sqstm1)) %>%
  dplyr::mutate(p62 = ifelse(p62 < -2, -2, p62)) %>%
  dplyr::mutate(p62 = ifelse(p62 > 7.5, 7.5, p62)) %>%
  ggplot(aes(y = p62, x = sqstm1, color = cancer_types)) +
  # coord_fixed(ratio = 1)
  geom_point() +
  scale_color_manual(
    values = p62_sqstm1_coef_label$color, 
    limits = p62_sqstm1_coef_label$cancer_types,
    labels = p62_sqstm1_coef_label$label
  ) +
  labs(
    x = "mRNA expression (log2)",
    y = "Protein expression",
    title = glue::glue("p62 protein vs. mRNA, R = {signif(p62_sqstm1_total_corr[1], digits = 2)}, p = {p62_sqstm1_total_corr[2]}")
  ) +
  theme(
    # for panel
    panel.background = element_rect(fill = "white", colour = "black", size = 1.3),
    panel.grid = element_blank(),
    
    # for legend
    legend.title = element_blank(),
    
    legend.background = element_rect(fill = "transparent", colour = "transparent"),
    legend.margin = margin(t = 0, r = 0, b = 0, l = 0),
    
    legend.key = element_rect(fill = "transparent", colour = "transparent"),
    legend.key.size = unit(x = 0.01, units = "npc"),
    legend.text = element_text(size = rel(0.7), colour = "black"),
    
    legend.position = c(0.235, 0.75)
  ) +
  guides(colour = guide_legend(ncol = 2)) ->
  plot_p62_sqstm1_expr

ggsave(
  filename = "02-p62-protein-vs-mrna.pdf", 
  plot = plot_p62_sqstm1_expr, 
  device = "pdf", 
  path = path_out,
  width = 7, height = 6.42
)

#
#


# sqstm1 correlation with atg5 and atg7 -----------------------------------
# atg5 & atg7 becn1 lc3
gene_set <- c("ATG5", "ATG7", "BECN1", "MAP1LC3A", "MAP1LC3B", "MAP1LC3C", "GABARAPL1", "GABARAPL2", "GABARAPL3", "LAMP2", "SQSTM1")
mrna_expr %>% 
  dplyr::mutate(
    gene_expr = purrr::map(
      .x = expr,
      .f = function(.x) {
        .x %>% 
          dplyr::filter(symbol %in% gene_set) %>% 
          dplyr::select(-2) -> .xx
        # select tumor
        names(.xx)[-1] -> .barcode
        .barcode[substr(x = .barcode, start = 14, stop = 14) == 0] -> .tumor
        
        .xx %>% 
          dplyr::select(1, .tumor) %>% 
          tidyr::gather(key = "barcode", value = "expr", -symbol) %>% 
          tidyr::spread(key = symbol, value = expr) %>% 
          dplyr::filter_if(.predicate = is.numeric, .vars_predicate = dplyr::all_vars(. != 0)) 
      }
    )
  ) %>% 
  dplyr::select(-2) %>% 
  tidyr::unnest() %>% 
  dplyr::select(-2) -> 
  genes_mrna_expr

# total correlation
combn(x = gene_set, m = 2, simplify = F) -> gene_set_comb

gene_set_comb %>% 
  purrr::map(
    .f = function(.x, .d = genes_mrna_expr) {
      cor.test(x = .d[[.x[1]]], y = .d[[.x[2]]], method = "pearson") %>% 
        broom::tidy() %>% 
        dplyr::select(coef = estimate, pval = p.value) %>% 
        tibble::add_column(vs = paste0(.x, collapse = "#"), .before = 1)
    }
  ) %>% 
  dplyr::bind_rows() %>% 
  dplyr::arrange(-coef) ->
  gene_set_total_corr

gene_set_total_corr %>% 
  dplyr::filter(grepl(pattern = "SQSTM1", x = vs)) %>%
  dplyr::arrange(coef)

# each cancer types
genes_mrna_expr %>% 
  dplyr::group_by(cancer_types) %>% 
  dplyr::do(
    gene_set_comb %>% 
      purrr::map(
        .f = function(.x, .d = .data) {
          cor.test(x = .d[[.x[1]]], y = .d[[.x[2]]], method = "pearson") %>% 
            broom::tidy() %>% 
            dplyr::select(coef = estimate, pval = p.value) %>% 
            tibble::add_column(vs = paste0(.x, collapse = "#"), .before = 1)
        }
      ) %>% 
      dplyr::bind_rows()
  ) %>% 
  dplyr::ungroup() ->
  gene_set_corr

gene_set_corr %>% 
  dplyr::filter(grepl(pattern = "SQSTM1", x = vs)) %>%
  dplyr::mutate(vs = sub(pattern = "#SQSTM1", replacement = "", x = vs, fixed = TRUE)) %>% 
  dplyr::rename(gene = vs) ->
  gene_set_corr_sqstm1

gene_set_corr_sqstm1 %>% 
  dplyr::group_by(cancer_types) %>% 
  dplyr::summarise(s = sum(coef)) %>% 
  dplyr::arrange(s) %>% 
  dplyr::pull(cancer_types) -> 
  rank_cancer

gene_set_corr_sqstm1 %>% 
  dplyr::group_by(gene) %>% 
  dplyr::summarise(s = sum(coef)) %>% 
  dplyr::arrange(s) %>% 
  dplyr::pull(gene) ->
  rank_gene

rank_gene <- rev(gene_set[-length(gene_set)])

gene_set_corr_sqstm1 %>% 
  dplyr::mutate(coef = ifelse(coef > 0.3, 0.3, coef)) %>%
  dplyr::mutate(coef = ifelse(coef < -0.3, -0.3, coef)) %>% 
  ggplot(aes(x = cancer_types, y = gene, fill = coef)) +
  geom_tile() +
  coord_fixed(ratio = 1) +
  scale_fill_gradient2(
    breaks = round(seq(-0.3, 0.3,length.out = 5), digits = 2),
    labels = format(seq(-0.3, 0.3,length.out = 5), digits = 2),
    low = "#33cbfe",
    mid = "#000000",
    high = "#fdfe00"
  ) +
  scale_x_discrete(limits = rank_cancer) +
  scale_y_discrete(limits = rank_gene) +
  labs(
    x = "",
    y = "",
    title = "p62 mRNA correlates with autophagy core genes."
  ) +
  theme(
    panel.background = element_blank(),
    panel.grid = element_blank(),
    
    # ticks
    axis.ticks = element_blank(),
    
    # axis text
    axis.text = element_text(size = 12, colour = "black"),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    
    # legend
    
    legend.position = "right"
  ) +
  guides(
    fill = guide_legend(
      # legend title
      title = "Pearsons' correlation (r)",
      title.position = "right",
      title.theme = element_text(angle = -90, size = 10),
      title.vjust = -0.3,
      title.hjust = 0.7,
      
      # legend label
      label.position = "left",
      label.theme = element_text(size = 14),
      
      # legend key
      keyheight = unit(x = 0.085, units = "npc"),
      
      reverse = TRUE
    )
  ) -> 
  plot_gene_set_corr_sqstm1

ggsave(
  filename = "03-p62-mrna-correlates-with-autophagy-flux-core-gene.pdf",
  plot = plot_gene_set_corr_sqstm1,
  device = "pdf",
  path = path_out,
  width = 10
)

#
#

# p62 correlates with pathway score ---------------------------------------
# p62 rppa expression correlates with the rppa score. ----

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
        .corrcomb %>% 
          purrr::map(.f = function(.cn){
            cor.test(x = .d[[.cn[1]]], y = .d[[.cn[2]]], method = "pearson") %>% 
              broom::tidy() %>% 
              dplyr::select(coef = estimate, pval = p.value) %>% 
              tibble::add_column(vs = paste0(.cn, collapse = "#"), .before = 1) 
          }) %>% 
          dplyr::bind_rows() 
      }
    )
  ) %>% 
  dplyr::select(cancer_types, corr) %>% 
  tidyr::unnest() -> 
  rppa_pathway_corr

# pathway name 
pathway_name <- c(
  "DNADamage" = "DNA Damage Response", 
  "CellCycle" = "Cell Cycle", 
  "TSCmTOR" = "TSC/mTOR",
  "PI3KAKT" = "PI3K/AKT",
  "RASMAPK" = "RAS/MAPK"
)

rppa_pathway_corr %>% 
  dplyr::filter(grepl(pattern = "#p62", x = vs)) %>% 
  dplyr::mutate(vs = sub(pattern = "#p62", replacement = "", x = vs)) %>% 
  dplyr::mutate(vs = plyr::revalue(x = vs, replace = pathway_name, warn_missing = F)) ->
  p62_pathway_corr

p62_pathway_corr %>% 
  dplyr::group_by(cancer_types) %>% 
  dplyr::summarise(m = sum(coef)) %>% 
  dplyr::arrange(m) %>% 
  dplyr::pull(cancer_types) ->
  rank_cancer

rank_cancer_pathway <- rank_cancer

p62_pathway_corr %>% 
  dplyr::group_by(vs) %>% 
  dplyr::summarise(m = sum(coef)) %>% 
  dplyr::arrange(-m) %>% 
  dplyr::pull(vs) ->
  rank_pathway

p62_pathway_corr %>% 
  dplyr::mutate(coef = ifelse(coef > 0.5, 0.5, coef)) %>%
  dplyr::mutate(coef = ifelse(coef < -0.5, -0.5, coef)) %>% 
  ggplot(aes(x = cancer_types, y = vs, fill = coef)) +
  geom_tile() +
  coord_fixed(ratio = 1) +
  scale_fill_gradient2(
    breaks = round(seq(-0.5, 0.5,length.out = 5), digits = 2),
    labels = format(seq(-0.5, 0.5,length.out = 5), digits = 2),
    low = "#00fefe",
    mid = "#000000",
    high = "#fe0000"
  ) +
  scale_x_discrete(limits = rank_cancer) +
  scale_y_discrete(limits = rank_pathway) +
  labs(
    x = "",
    y = "",
    title = "p62 protein expresion correlates with rppa score"
  ) +
  theme(
    panel.background = element_blank(),
    panel.grid = element_blank(),
    
    # ticks
    axis.ticks = element_blank(),
    
    # axis text
    axis.text = element_text(size = 12, colour = "black"),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    
    # legend
    
    legend.position = "right"
  ) +
  guides(
    fill = guide_legend(
      # legend title
      title = "Pearsons' correlation (r)",
      title.position = "right",
      title.theme = element_text(angle = -90, size = 11),
      title.vjust = -0.3,
      title.hjust = 0.7,
      
      # legend label
      label.position = "left",
      label.theme = element_text(size = 12),
      
      # legend key
      keyheight = unit(x = 0.08, units = "npc"),
      
      reverse = TRUE
    )
  ) -> 
  plot_p62_pathway_corr

ggsave(
  filename = "04-p62-correlates-with-rppa-pathway.pdf",
  plot = plot_p62_pathway_corr,
  device = "pdf",
  path = path_out,
  width = 10
)
##
#

# protein in the pathway -----
rppa_expr %>% 
  dplyr::mutate(
    rppa = purrr::map2(
      .x = expr,
      .y = cancer_types,
      .f = function(.x, .y) {
        print(.y)
        pathway_gene_list %>% 
          dplyr::select(-3) %>% 
          tidyr::unnest() %>% 
          dplyr::mutate(protein = gsub("-R-C|-R-V|-M-C|-M-V|-R-E|-G-C|-R-E|-M-E", "", protein)) %>% 
          dplyr::mutate(protein = ifelse(startsWith(x = protein, prefix = "4"), glue::glue("X{protein}"), protein)) %>% 
          dplyr::mutate(protein = ifelse(startsWith(x = protein, prefix = "5"), glue::glue("X{protein}"), protein)) %>% 
          dplyr::mutate(pat = gsub(pattern = "-|_", replacement = "", x = protein) %>% tolower()) %>% 
          # c("retpy905", "caspase9", "parp1", "parpab3")
          dplyr::filter(!pat %in% c("retpy905", "caspase9", "parp1", "parpab3")) %>% 
          dplyr::select(pathway, pat) ->
          .pathway_gene_list
        
        .x %>% tibble::add_column(pat = gsub(pattern = "-|_", replacement = "", x = .x$protein) %>% tolower(), .before = 3) -> .xx
        .xx %>% dplyr::filter(protein == "P62LCKLIGAND") -> .p62
        # setdiff(.pathway_gene_list$pat, .xx$pat)
        
        .xx %>% dplyr::filter(pat %in% .pathway_gene_list$pat) -> .d
        
        .d %>% 
          dplyr::bind_rows(.p62) %>% 
          dplyr::select(-symbol, -protein) %>% 
          tidyr::gather(key = "barcode", value = "rppa", -pat) %>% 
          tidyr::spread(key = pat, value = rppa) -> 
          .dd
        
        .d$pat %>% sort() %>% unique() %>% purrr::map(.f = function(.x) {c(.x, "p62lckligand")}) -> .corrcomb
        
        .corrcomb %>% 
          purrr::map(.f = function(.cn){
            cor.test(x = .dd[[.cn[1]]], y = .dd[[.cn[2]]], method = "pearson") %>% 
              broom::tidy() %>% 
              dplyr::select(coef = estimate, pval = p.value) %>% 
              tibble::add_column(pat = paste0(.cn, collapse = "#"), .before = 1)
        }) %>% 
          dplyr::bind_rows() %>% 
          dplyr::mutate(pat = sub(pattern = "#p62lckligand", replacement = "", x = pat)) ->
          .corr
        
        .pathway_gene_list %>% 
          dplyr::inner_join(.xx, by = "pat") %>% 
          dplyr::select(pathway, pat, protein) %>% 
          dplyr::inner_join(.corr, by = "pat") %>% 
          dplyr::select(-pat)
      }
    )
  ) ->
  p62_protein_pathway_corr

p62_protein_pathway_corr %>% 
  dplyr::select(-expr) %>% 
  tidyr::unnest() ->
  p62_protein_pathway_corr_p


the_pathway <- c("PI3K/AKT", "RAS/MAPK", "RTK", "TSC/mTOR", "EMT")

the_pathway %>% 
  purrr::map(
    .f = function(.x) {
      
      p62_protein_pathway_corr_p %>% 
        dplyr::filter(pathway == .x) ->
        .path
      
      # .path %>% 
      #   dplyr::group_by(cancer_types) %>% 
      #   dplyr::summarise(m = sum(coef)) %>% 
      #   dplyr::arrange(m) %>% 
      #   dplyr::pull(cancer_types) ->
      #   rank_cancer
      
      .path %>% 
        dplyr::group_by(protein) %>% 
        dplyr::summarise(m = sum(coef)) %>% 
        dplyr::arrange(m) %>% 
        dplyr::pull(protein) ->
        rank_protein
      
      .path %>% 
        dplyr::mutate(
          coef = dplyr::case_when(
            coef > 0.5 ~ 0.5,
            coef < -0.5 ~ -0.5,
            TRUE ~ coef
          )
        ) %>% 
        ggplot(aes(x = cancer_types, y = protein, fill = coef)) +
        geom_tile() +
        coord_fixed(ratio = 1) +
        scale_fill_gradient2(
          breaks = round(seq(-0.5, 0.5,length.out = 5), digits = 2),
          labels = format(seq(-0.5, 0.5,length.out = 5), digits = 2),
          low = "#00fefe", # 1cb5e2, 33cbfe, 00fefe
          mid = "#000000",
          high = "#fe0000" # fdfe00, fe00fe, fe0000
        ) +
        scale_x_discrete(limits = rank_cancer_pathway) +
        scale_y_discrete(limits = rank_protein, labels = gsub(pattern = "_", replacement = " ", x = rank_protein)) +
        labs(
          x = "",
          y = "",
          title = glue::glue("p62 correlates with {.x} pathway proteins")
        ) +
        theme(
          panel.background = element_blank(),
          panel.grid = element_blank(),
          
          # ticks
          axis.ticks = element_blank(),
          
          # axis text
          axis.text = element_text(size = 12, colour = "black"),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
          
          # legend
          
          legend.position = "right"
        ) +
        guides(
          fill = guide_legend(
            # legend title
            title = "Pearsons' correlation (r)",
            title.position = "right",
            title.theme = element_text(angle = -90, size = 11),
            title.vjust = -0.3,
            title.hjust = 0.7,
            
            # legend label
            label.position = "left",
            label.theme = element_text(size = 12),
            
            # legend key
            keyheight = unit(x = 0.08, units = "npc"),
            
            reverse = TRUE
          )
        ) -> 
        .plot
      .filename <- glue::glue("{sub(pattern = '/', replacement = '-', x = .x)}-protein-p62-correlation.pdf")
      .path_plot_pathway <- file.path(path_out, "02-p62-correlates-with-pathway-protein")
      if (!dir.exists(.path_plot_pathway)) dir.create(path = .path_plot_pathway)
      ggsave(
        filename = .filename,
        plot = .plot,
        device = "pdf",
        path = .path_plot_pathway,
        width = 10
      )
    }
  )


