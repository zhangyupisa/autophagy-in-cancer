library(magrittr)
library(ggplot2)

expr_path <- "/home/cliu18/liucj/projects/6.autophagy/02_autophagy_expr/"
stage_path <- file.path(expr_path, "03_b_stage")
tcga_path <- "/home/cliu18/liucj/projects/6.autophagy/TCGA"
expr_path <- file.path(expr_path, "03_a_gene_expr")

clinical_stage <- 
  readr::read_rds(path = file.path(tcga_path,"pancan34_clinical_stage.rds.gz")) %>% 
  dplyr::filter(n >= 40) %>% 
  dplyr::select(-n)

gene_list <- readr::read_rds(file.path(expr_path, "rds_03_a_atg_lys_gene_list.rds.gz"))
gene_list_expr <- readr::read_rds(path = file.path(expr_path, ".rds_03_a_gene_list_expr.rds.gz"))


# merge clinical and expr
expr_stage <- 
  gene_list_expr %>%
  dplyr::inner_join(clinical_stage, by = "cancer_types")

fun_barcode <- function(.b){
  stringr::str_sub(
    string = .b,
    start = 1,
    end = 12
  )
} #get short barcode from long barcode
fun_tn_type <- function(.b){
  type <- .b %>% 
    stringr::str_split(pattern = "-", simplify = T) %>% 
    .[, 4] %>% 
    stringr::str_sub(1, 2)
} # get tumor and normal info
fun_expr_stage_merge <- function(filter_expr, stage){
  # merge clinical and expr data
  filter_expr %>% 
    dplyr::select(-entrez_id) %>% 
    tidyr::gather(key = barcode, value = expr, -symbol) %>% 
    dplyr::mutate(type = fun_tn_type(barcode)) %>% 
    dplyr::filter(type == "01") %>% 
    dplyr::mutate(barcode = fun_barcode(barcode)) %>% 
    dplyr::select(symbol, barcode, expr) -> expr_clean
  expr_clean %>% dplyr::inner_join(stage, by = "barcode") -> expr_stage_ready
} # merge stage and expr to clean
fn_get_order <- function(.d){
  .d %>% 
    dplyr::group_by(stage ) %>% 
    dplyr::summarise(me = mean(expr)) %>% 
    .$me %>% rank() -> .d_m
  
  if(identical(.d_m, c(1,2,3,4))){
    return(1)
  } else if(identical(.d_m, c(4,3,2,1))){
    return(2)
  } else{
    return(3)
  }
}
fun_stage_test <- function(expr_stage_ready){
  expr_stage_ready %>% 
    tidyr::drop_na(expr) %>%
    dplyr::group_by(symbol, stage) %>% 
    dplyr::mutate(l = n() > 10) -> tl
  
  if(! all(tl$l)){
    return(tibble::tibble())
  } else{
    expr_stage_ready %>% 
      tidyr::drop_na(expr) %>%
      dplyr::group_by(symbol) %>% 
      dplyr::do(broom::tidy(oneway.test(log2(expr + 1) ~ stage, data = .))) %>%
      dplyr::ungroup() %>% 
      dplyr::mutate(fdr = p.adjust(p.value, method = "fdr")) %>% 
      dplyr::select(symbol, p.value, fdr) -> diff_pval
    
    expr_stage_ready %>% 
      tidyr::drop_na(expr) %>%
      tidyr::nest(-symbol) %>% 
      dplyr::mutate(order = purrr::map_dbl(data, .f = fn_get_order)) %>% 
      dplyr::select(- data) -> symbol_order

    diff_pval %>% 
      dplyr::inner_join(symbol_order, by = "symbol")
  }
} # stage_test


# expr_stage %>%
#   dplyr::filter(cancer_types == "BRCA") %>% 
#   dplyr::mutate(merged_clean = purrr::map2(filter_expr, stage, fun_expr_stage_merge)) %>%
#   dplyr::select(-filter_expr, -stage) %>%
#   dplyr::mutate(diff_pval = purrr::map(merged_clean, fun_stage_test)) %>%
#   tidyr::unnest(diff_pval, .drop = F) -> te

cl <- parallel::detectCores()
cluster <- multidplyr::create_cluster(floor(cl * 5 / 6))
expr_stage %>%
  multidplyr::partition(cluster = cluster) %>%
  multidplyr::cluster_library("magrittr") %>%
  multidplyr::cluster_assign_value("fun_barcode", fun_barcode)  %>%
  multidplyr::cluster_assign_value("fun_tn_type", fun_tn_type) %>% 
  multidplyr::cluster_assign_value("fun_expr_stage_merge", fun_expr_stage_merge) %>% 
  multidplyr::cluster_assign_value("fun_stage_test", fun_stage_test) %>% 
  multidplyr::cluster_assign_value("fn_get_order", fn_get_order) %>% 
  dplyr::mutate(merged_clean = purrr::map2(filter_expr, stage, fun_expr_stage_merge)) %>% 
  dplyr::select(-filter_expr, -stage) %>% 
  dplyr::mutate(diff_pval = purrr::map(merged_clean, fun_stage_test)) %>% 
  dplyr::collect() %>%
  dplyr::as_tibble() %>%
  dplyr::ungroup() %>%
  dplyr::select(-PARTITION_ID) %>% 
  tidyr::unnest(diff_pval, .drop = F) -> expr_stage_sig_pval
on.exit(parallel::stopCluster(cluster))

expr_stage_sig_pval %>% 
  readr::write_rds(path = file.path(stage_path, ".rds_03_b_stage_gene_expr.rds.gz"), compress = "gz")
# expr_stage_sig_pval <- readr::read_rds(path = file.path(stage_path, ".rds_03_b_stage_gene_expr.rds.gz"))
#--------------------------------------------------------

fun_rank_cancer <- function(pattern){
  pattern %>% 
    dplyr::summarise_if(.predicate = is.numeric, dplyr::funs(sum(., na.rm = T))) %>%
    tidyr::gather(key = cancer_types, value = rank) %>%
    dplyr::arrange(dplyr::desc(rank))
} #get cancer rank
fun_rank_gene <- function(pattern){
  pattern %>% 
    dplyr::rowwise() %>%
    dplyr::do(
      symbol = .$symbol,
      rank =  unlist(.[-1], use.names = F) %>% sum(na.rm = T)
    ) %>%
    dplyr::ungroup() %>%
    tidyr::unnest() %>%
    dplyr::arrange(rank)
} # get gene rank


expr_stage_sig_pval %>% 
  dplyr::select(cancer_types, symbol) %>% 
  dplyr::mutate(n = 1) %>% 
  tidyr::spread(key = cancer_types, value = n) -> pattern

cancer_rank <- pattern %>% fun_rank_cancer()
gene_rank <- pattern %>% 
  fun_rank_gene() %>% 
  dplyr::left_join(gene_list, by = "symbol") %>% 
  dplyr::mutate(color = plyr::revalue(status, replace = c('a' = "#e41a1c", "l" = "#377eb8", "i" = "#4daf4a", "p" = "#984ea3"))) %>%
  dplyr::arrange(color, rank)

expr_stage_sig_pval %>% 
  ggplot(aes(x = cancer_types, y = symbol, color = cancer_types)) +
  geom_point(aes(size = -log10(p.value))) +
  scale_x_discrete(limit = cancer_rank$cancer_types) +
  scale_y_discrete(limit = gene_rank$symbol) +
  scale_size_continuous(
    limit = c(-log10(0.05), 15),
    range = c(1, 6),
    breaks = c(-log10(0.05), 5, 10, 15),
    labels = c("0.05", latex2exp::TeX("$10^{-5}$"), latex2exp::TeX("$10^{-10}$"), latex2exp::TeX("$< 10^{-15}$")),
    name = "P-value"
  ) +
  theme(
    panel.background = element_rect(colour = "black", fill = "white"),
    panel.grid = element_line(colour = "grey", linetype = "dashed"),
    panel.grid.major = element_line(
      colour = "grey",
      linetype = "dashed",
      size = 0.2
    ),
    
    axis.title = element_blank(),
    axis.ticks = element_line(color = "black"),
    axis.text.y = element_text(color = gene_rank$color),
    
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14),
    legend.key = element_rect(fill = "white", colour = "black")
  ) ->p
ggsave(
  filename = "fig_03_b_stage_sig_genes.pdf",
  plot = p,
  device = "pdf",
  width = 8,
  height = 22,
  path = stage_path
)
readr::write_rds(
  p,
  path = file.path(stage_path, ".fig_03_b_stage_sig_genes.pdf.rds.gz"),
  compress = "gz"
)
#-------------------------------------
  
fun_draw_boxplot <- function(cancer_types, merged_clean, symbol, p.value, fdr){
  # print(cancer_types)
  p_val <- signif(-log10(p.value), digits = 3)
  gene <- symbol
  fig_name <- paste(cancer_types, gene, p_val, "pdf", sep = ".")
  comp_list <- list(c("Stage I", "Stage II"), c("Stage II", "Stage III"), c("Stage III", "Stage IV"))
  merged_clean %>% 
    dplyr::filter(symbol == gene) %>% 
    dplyr::mutate(expr = log2(expr)) %>% 
    dplyr::arrange(stage) %>% 
    ggpubr::ggboxplot(x = "stage", y = "expr",  color = "stage", pallete = "jco"  ) +
    ggpubr::stat_compare_means(comparisons = comp_list, method = "t.test") + 
    ggpubr::stat_compare_means(method = "anova", label.y = 14) +
    labs(x  = "", y = "Expression (log2 RSEM)", title = paste(gene, "expression stage change in", cancer_types)) +
    ggthemes::scale_color_gdocs() -> p
    ggsave(filename = fig_name, plot = p, path = file.path(stage_path, "boxplot"), width = 6, height = 6,  device = "pdf")
}
fun_draw_boxplot_filter <- function(cancer_types, merged_clean, symbol, p.value, fdr){
  p_val <- signif(-log10(p.value), digits = 3)
  gene <- symbol
  fig_name <- paste(cancer_types, gene, p_val, "pdf", sep = ".")
  comp_list <- list(c("Stage I", "Stage II"), c("Stage II", "Stage III"), c("Stage III", "Stage IV"))
  
  d <- merged_clean %>% 
    dplyr::filter(symbol == gene) %>% 
    dplyr::mutate( expr = log2(expr)) %>% 
    dplyr::arrange(stage)
    
  p_list <- ggpubr::compare_means(expr ~ stage, data = d, method = "t.test")
  
  if(p_list %>% .$p %>% .[3] < 0.05){
    d %>% 
      ggpubr::ggboxplot(x = "stage", y = "expr",  color = "stage", pallete = "jco"  ) +
      ggpubr::stat_compare_means(comparisons = comp_list, method = "t.test") + 
      ggpubr::stat_compare_means(method = "anova", label.y = 14) +
      labs(x  = "", y = "Expression (log2 RSEM)", title = paste(gene, "expression subtype change in", cancer_types)) +
      ggthemes::scale_color_gdocs() -> p
      ggsave(filename = fig_name, plot = p, path = file.path(stage_path, "boxplot_stageI_vs_IV"), width = 6, height = 6,  device = "pdf")
  } else{
    print(fig_name)
  }
}


expr_stage_sig_pval %>% purrr::pwalk(.f = fun_draw_boxplot)
# expr_stage_sig_pval %>% purrr::pwalk(.f = fun_draw_boxplot_filter)


save.image(file = file.path(stage_path, ".rda_03_b_stage_gene_expr.rda"))
load(file = file.path(stage_path, ".rda_03_b_stage_gene_expr.rda"))
