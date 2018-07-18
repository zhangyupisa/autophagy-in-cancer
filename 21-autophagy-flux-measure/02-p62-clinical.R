
# Load library ------------------------------------------------------------

library(magrittr)
library(ggplot2)

# Path --------------------------------------------------------------------

path_project <- "/home/liucj/data/project/06-autophagy"
path_data <- file.path(path_project, "TCGA")
path_out <- file.path(path_project, "21-autophagy-flux")


# p62 survival ------------------------------------------------------------

# load data ---------------------------------------------------------------
# pancan33-rppa-expr-v4-l4.rds.gz
rppa_expr <- readr::read_rds(path = file.path(path_data, "pancan33-rppa-expr-v4-l4.rds.gz")) %>% 
  dplyr::filter(cancer_types != "PANCAN19")

# TCGA_pancan_cancer_cell_survival_time.rds.gz
clinical <- readr::read_rds(path = file.path(path_data, "TCGA_pancan_cancer_cell_survival_time.rds.gz")) %>% 
  dplyr::mutate(data = purrr::map(
    .x = data, 
    .f = function(.x) .x %>% dplyr::select("barcode" = "bcr_patient_barcode", "status" = "PFI.1", "time" = "PFI.time.1"))
  ) %>% 
  dplyr::rename(cancer_types = type)
tcga_cell_clinical <- file.path(path_data, "tcga-cell-survival-lite.rds.gz")
if (!file.exists(tcga_cell_clinical)) readr::write_rds(clinical, path = tcga_cell_clinical, compress = "gz")

# p62 rppa expresion ------------------------------------------------------

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

# merge p62 expr with survival data
clinical %>% 
  dplyr::inner_join(p62_rppa_expr, by = "cancer_types") %>% 
  dplyr::mutate(
    survival = purrr::map2(
      .x = data,
      .y = p62,
      .f = function(.x, .y) {
        .y %>% 
          tidyr::gather(key = "barcode", value = "p62", -protein) %>% 
          dplyr::select(-protein) %>% 
          dplyr::mutate(barcode = substr(x = barcode, start = 1, stop = 12)) %>% 
          dplyr::group_by(barcode) %>% 
          dplyr::summarise(p62 = mean(p62)) %>% 
          dplyr::inner_join(.x, by = "barcode") %>% 
          dplyr::filter(!is.na(time), time > 0, !is.na(status)) %>% 
          dplyr::mutate(time = ifelse(time / 30 > 60, 60, time / 30)) %>%
          dplyr::mutate(status = ifelse(time >= 60 & status == 1, 0, status)) -> 
          .d
      }
    )
  ) %>% 
  dplyr::select(-p62, -data) ->
  p62_survival

p62_survival %>% 
  dplyr::mutate(
    hazard = purrr::map2(
      .x = cancer_types,
      .y = survival,
      .f = function(.x, .y) {
        .y %>% 
          dplyr::mutate(
            group = dplyr::case_when(
              p62 >= quantile(p62, 0.5) ~ "L",
              p62 < quantile(p62, 0.5) ~ "H",
              TRUE ~ "M"
            )
          ) -> 
          .d
        
        survival::coxph(survival::Surv(time, status) ~ p62, data = .d) -> .cox
        summary(.cox) -> .z
        survival::survdiff(survival::Surv(time, status) ~ group, data = .d) -> .diff
        .kmp <- 1 - pchisq(.diff$chisq, df = length(levels(as.factor(.d$group))) - 1)
        
        tibble::tibble(
          n = .z$n,
          hr = .z$conf.int[1],
          hr_l = .z$conf.int[3],
          hr_h = .z$conf.int[4],
          coxp = .z$waldtest[3],
          kmp = .kmp
          )
      }
    )
  ) %>% 
  dplyr::select(-survival) %>% 
  tidyr::unnest() -> 
  p62_hazard_ratio


p62_hazard_ratio %>% 
  dplyr::mutate(coxp = -log10(coxp)) %>% 
  dplyr::mutate(cancer_types = glue::glue("{cancer_types} ({n})")) %>% 
  dplyr::mutate(cancer_types = reorder(cancer_types, hr, sort)) -> 
  plot_ready

plot_ready %>% 
  ggplot(aes(y = hr, x = cancer_types)) +
  geom_point(aes(size = coxp)) +
  geom_hline(aes(yintercept = 1), colour = "red") +
  # geom_hline(aes(yintercept = 1.1), colour = "green", linetype = "dashed") +
  # geom_hline(aes(yintercept = 0.9), colour = "green", linetype = "dashed") +
  scale_size(name = "p-value") +
  coord_flip() +
  labs(y = "Hazard Ratio", x = "Cancer Types") +
  theme(
    panel.background = element_blank(),
    panel.border = element_rect(colour = "black", fill = "transparent"),
    text = element_text(family = "serif")
  ) -> 
  plot_p62_hazard_ratio

ggsave(
  filename = "05-p62-hazard-ratio.pdf",
  plot = plot_p62_hazard_ratio,
  device = "pdf",
  path = path_out
)

knitr::kable(plot_ready %>% 
               dplyr::rename(`hazard ratio` = hr) %>% 
               dplyr::mutate(coxp = 10^-coxp) %>% 
               dplyr::arrange(-`hazard ratio`) %>% 
               dplyr::mutate(hr_l = signif(hr_l, 3), hr_h = signif(hr_h,3)) %>% 
               tidyr::unite(`95%CI`, hr_l, hr_h, sep = "-"))


# p62 survival plot -------------------------------------------------------
p62_survival %>% 
  dplyr::mutate(
    hazard = purrr::map2(
      .x = cancer_types,
      .y = survival,
      .f = function(.x, .y) {
        .y %>% 
          dplyr::mutate(
            group = dplyr::case_when(
              p62 >= quantile(p62, 0.5) ~ "L",
              p62 < quantile(p62, 0.5) ~ "H",
              TRUE ~ "M"
            )
          ) -> 
          .d
        
        survival::coxph(survival::Surv(time, status) ~ p62, data = .d) %>% 
          broom::tidy() %>% 
          dplyr::mutate(hazard_ratio = exp(estimate)) %>% 
          dplyr::select(hazard_ratio, coxp = p.value) -> 
          .hazard_coxp
        
        survival::survdiff(survival::Surv(time, status) ~ group, data = .d) -> .diff
        .kmp <- 1 - pchisq(.diff$chisq, df = length(levels(as.factor(.d$group))) - 1)
        
        survival::survfit(survival::Surv(time, status) ~ group, data = .d , na.action = na.exclude) -> .fit
        
        CPCOLS <- c("#CD0000", "#000080")
        survminer::ggsurvplot(
          .fit, data = .d, 
          # pval = TRUE, pval.method = T,
          xlab = "Survival in months", ylab = 'Probability of survival',
          palette = CPCOLS, ggtheme = theme_bw()
          ) -> 
          .plot
        
        .label <- glue::glue("{.x} 
                             Log-rank p = {signif(.kmp, 2)}")
        .plot$plot +
          annotate(
            "text", 
            x = ifelse(.x == "GBM", 40, 10), 
            y = ifelse(.x == "GBM", 0.8, 0.2), 
            label = .label
            ) + 
          theme(
            legend.position = "",
            panel.grid = element_blank()
            ) ->
          .plot_label
        
        tibble::tibble(
          hazard_ratio = .hazard_coxp$hazard_ratio,
          coxp = .hazard_coxp$coxp,
          kmp = .kmp,
          plot = list(.plot_label)
        )
      }
    )
  ) %>% 
  dplyr::select(-survival) %>% 
  tidyr::unnest() ->
  p62_coxp_kmp_plot

gridExtra::arrangeGrob(
  grobs = p62_coxp_kmp_plot %>% dplyr::filter(kmp < 0.05) %>% .$plot, 
  nrow = 2
  ) %>% 
  ggsave(
    filename = "06-p62-cancer-type-survival.pdf",
    plot = .,
    device = "pdf", 
    width = 17,
    height = 8,
    path = path_out
  )

#
#
# end  --------------------------------------------------------------------




