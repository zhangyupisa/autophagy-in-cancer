
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
            x = ifelse(.x == "GBM", 45, 15), 
            y = ifelse(.x == "GBM", 0.75, 0.2), 
            label = .label
            ) + 
          theme(
            text = element_text(family = "serif"),
            legend.position = "none",
            # legend.title = element_blank(),
            # legend.background = element_rect(colour = "black", fill = "transparent"),
            # legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"),
            # legend.key = element_rect(colour = "black", fill = "transparent"),
            axis.title = element_blank(),
            
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
  nrow = 2, bottom = "Survival in months", left = "Probability of survival"
  ) %>% 
  ggsave(
    filename = "06-p62-cancer-type-survival.pdf",
    plot = .,
    device = "pdf",
    width = 17,
    height = 8,
    path = path_out
  )


# stage info for the p62 --------------------------------------------------


# load stage data ---------------------------------------------------------
# pancan34-clinical-stage.rds.gz
stage <- readr::read_rds(path = file.path(path_data, "pancan34-clinical-stage.rds.gz")) %>% 
  dplyr::filter(n >= 40) %>% 
  dplyr::select(-n)
# each stage has 10 samples at least.

p62_rppa_expr %>% 
  dplyr::inner_join(stage, by = "cancer_types") %>% 
  dplyr::mutate(
    p62_gender = purrr::pmap(
      .l = .,
      .f = function(cancer_types, p62, stage) {
        .x <- p62
        .y <- stage
        .z <- cancer_types
        
        .x %>% 
          tidyr::gather(key = "barcode", value = "p62", -protein) %>% 
          dplyr::select(-protein) %>% 
          dplyr::mutate(barcode = substr(x = barcode, start = 1, stop = 12)) %>% 
          dplyr::group_by(barcode) %>% 
          dplyr::summarise(p62 = mean(p62)) %>% 
          dplyr::inner_join(.y, by = "barcode") %>% 
          dplyr::select(barcode, p62, stage) %>% 
          tidyr::drop_na() %>% 
          dplyr::mutate(p62 = ifelse(p62 > 4, 4, p62)) ->
          .d
        
        .comp_list <- list(c("Stage I", "Stage II"), c("Stage II", "Stage III"), c("Stage III", "Stage IV"))
        
        .d %>% 
          dplyr::arrange(stage) %>% 
          ggpubr::ggboxplot(x = "stage", y = "p62",  color = "stage", pallete = "jco", add = "jitter") +
          ggpubr::stat_compare_means(comparisons = .comp_list, method = "wilcox.test") + 
          ggpubr::stat_compare_means(
            method = "kruskal.test", 
            label.x = 1.6,
            label.y = max(.d$p62) + 2.5, 
            aes(label = paste0(.z, ", Kruskal-Wallis, p = ", ..p.format..))
          ) +
          scale_color_brewer(palette = "Set1") +
          theme(
            text = element_text(family = "serif"),
            
            panel.background = element_rect(fill = "white", colour = "black", size = 1),
            axis.title = element_blank(),
            axis.line = element_blank(),
            axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
            
            legend.position = "none"
          ) ->
          .plot
        
        kruskal.test(p62 ~ as.factor(stage), data = .d) %>% 
          broom::tidy() %>% 
          dplyr::pull(p.value) ->
          .pval
        
        tibble::tibble(
          pval = .pval,
          plot = list(.plot)
        )
      }
    )
  ) %>% 
  dplyr::select(1, 4) %>% 
  tidyr::unnest() ->
  p62_stage

gridExtra::arrangeGrob(
  grobs = p62_stage %>% dplyr::filter(pval < 0.05) %>% .$plot, 
  nrow = 2, bottom = "Stage", left = "p62 protein expression"
  ) %>% 
  ggsave(
    filename = "07-p62-cancer-type-stage.pdf",
    plot = .,
    device = "pdf",
    width = 7,
    height = 7,
    path = path_out
  )


# gender ------------------------------------------------------------------
# pancan34-clinical.rds.gz
clinical_gender <- readr::read_rds(path = file.path(path_data, "pancan34-clinical.rds.gz"))
clinical_gender %>% 
  dplyr::mutate(
    gender = purrr::map(
      .x = clinical,
      .f = function(.x) {
        .x %>% 
          dplyr::select(barcode, gender) %>% 
          tidyr::drop_na(gender) %>% 
          dplyr::distinct() 
      }
    )
  ) %>% 
  dplyr::select(-clinical) %>% 
  dplyr::mutate(
    n = purrr::map(
      .x = gender,
      .f = function(.x) {
        .x$gender %>% table() -> .t
        tibble::tibble(
          F_n = ifelse(is.na(.t['FEMALE']), 0, .t['FEMALE']),
          M_n = ifelse(is.na(.t['MALE']), 0, .t['MALE'])
        )
      }
    )
  ) %>% 
  tidyr::unnest(n) ->
  gender

if (!file.exists(file.path(path_data, "pancan34-clinical-gender.rds.gz"))) 
  readr::write_rds(x = gender, path = file.path(path_data, "pancan34-clinical-gender.rds.gz"), compress = "gz")

# load gender data
readr::read_rds(path = file.path(path_data, "pancan34-clinical-gender.rds.gz")) %>% 
  dplyr::filter_if(.predicate = is.double, .vars_predicate = dplyr::all_vars(. > 10)) %>% 
  dplyr::select(1, 2) ->
  gender

p62_rppa_expr %>% 
  dplyr::inner_join(gender, by = "cancer_types") %>% 
  dplyr::mutate(
    p62_gender = purrr::pmap(
      .l = .,
      .f = function(cancer_types, p62, gender) {
        .x <- p62
        .y <- gender
        .z <- cancer_types
        
        print(.z)
        .x %>% 
          tidyr::gather(key = "barcode", value = "p62", -protein) %>% 
          dplyr::select(-protein) %>% 
          dplyr::mutate(barcode = substr(x = barcode, start = 1, stop = 12)) %>% 
          dplyr::group_by(barcode) %>% 
          dplyr::summarise(p62 = mean(p62)) %>% 
          dplyr::inner_join(.y, by = "barcode") %>% 
          tidyr::drop_na() %>% 
          dplyr::mutate(gender = factor(x = gender, levels = c("FEMALE", "MALE"))) ->
          .d
        
        .d$gender %>% table() %>% as.numeric() -> .dg
        
        if (gtools::invalid(.dg) || any(.dg < c(10, 10))) {
          return(
            tibble::tibble(
              f_m = NA,
              pval = NA,
              plot = list(NA)
            )
          )
          }
        
        t.test(p62 ~ gender, data = .d) %>% 
          broom::tidy() %>% 
          dplyr::select(f_m = estimate, f = estimate1, m = estimate2, pval = p.value ) ->
          .t
        
        .d %>% 
          ggpubr::ggboxplot(x = "gender", y = "p62",  color = "gender", pallete = "jco", add = "jitter") +
          ggpubr::stat_compare_means(
            method = "t.test", 
            label.x = 1,
            label.y = max(.d$p62) + 2.5, 
            aes(label = paste0(.z, ", T-test, p = ", ..p.format..))
          ) +
          scale_color_brewer(palette = "Set1") +
          theme(
            text = element_text(family = "serif"),
            
            panel.background = element_rect(fill = "white", colour = "black", size = 1),
            axis.title = element_blank(),
            axis.line = element_blank(),
            
            legend.position = "none"
          ) ->
          .plot
        
        tibble::tibble(
          f_m = .t$f_m,
          pval = .t$pval,
          plot = list(.plot)
        )
      }
    )
  ) %>% 
  dplyr::select(cancer_types, p62_gender) %>% 
  tidyr::unnest() ->
  p62_gender


gridExtra::arrangeGrob(
  grobs = p62_gender %>% dplyr::filter(pval < 0.05) %>% .$plot, 
  nrow = 2, bottom = "Gender", left = "p62 protein expression"
  ) %>% 
  ggsave(
    filename = "08-p62-cancer-type-gender.pdf",
    plot = .,
    device = "pdf",
    width = 7,
    height = 7,
    path = path_out
  )



# subtypes ----------------------------------------------------------------
# pancan34-clinical-subtype.rds.gz
subtype <- readr::read_rds(path = file.path(path_data, "pancan34-clinical-subtype.rds.gz")) %>% 
  dplyr::select(-n)


p62_rppa_expr %>% 
  dplyr::inner_join(subtype, by = "cancer_types") %>% 
  dplyr::mutate(
    p62_subtype = purrr::pmap(
      .l = .,
      .f = function(cancer_types, p62, subtype) {
        .x <- p62
        .y <- subtype
        .z <- cancer_types
        print(.z)
        
        .x %>% 
          tidyr::gather(key = "barcode", value = "p62", -protein) %>% 
          dplyr::select(-protein) %>% 
          dplyr::mutate(barcode = substr(x = barcode, start = 1, stop = 12)) %>% 
          dplyr::group_by(barcode) %>% 
          dplyr::summarise(p62 = mean(p62)) %>% 
          dplyr::inner_join(.y, by = "barcode") %>% 
          dplyr::select(barcode, p62, subtype) %>% 
          tidyr::drop_na()  ->
          .d
        
        .d %>% 
          dplyr::group_by(subtype) %>% 
          dplyr::mutate(l = n() >= 10) -> 
          .ds
        
        .ds$subtype %>% unique() %>% length() -> .ns
        
        .ds$subtype %>% unique() %>% sort() -> .lev
        
        if (!all(.ds$l) || .ns < 2) return(tibble::tibble(pval = NA, plot = list(NA)))
        
        # .ns == 2 t.test, .ns > 2 anova
        if (.ns == 2) {
          t.test(p62 ~ subtype, .d) %>% 
            broom::tidy() %>% 
            dplyr::pull(p.value) ->
            .pval
          
          .d %>% 
            ggpubr::ggboxplot(x = "subtype", y = "p62",  color = "subtype", pallete = "jco", add = "jitter") +
            ggpubr::stat_compare_means(
              method = "t.test", 
              label.x = 1,
              label.y = max(.d$p62) + 2.5, 
              aes(label = paste0(.z, ", T-test, p = ", ..p.format..))
            ) +
            scale_color_brewer(palette = "Set1") +
            theme(
              text = element_text(family = "serif"),
              
              panel.background = element_rect(fill = "white", colour = "black", size = 1),
              axis.title = element_blank(),
              axis.line = element_blank(),
              
              legend.position = "none"
            ) ->
            .plot
          
          tibble::tibble(
            pval = .pval,
            plot = list(.plot)
          )
        } else {
          oneway.test(p62 ~ subtype, .d) %>% 
            broom::tidy() %>% 
            dplyr::pull(p.value) -> 
            .pval
          
          .comp_list <- .lev %>% combn(m = 2, simplify = F)
          .d %>% dplyr::mutate(p62 = ifelse(p62 > 4, 4, p62)) -> .d
          .d %>% 
            dplyr::mutate(subtype = factor(subtype, levels = .lev)) %>% 
            ggpubr::ggboxplot(x = "subtype", y = "p62",  color = "subtype", pallete = "jco", add = "jitter") +
            ggpubr::stat_compare_means(comparisons = .comp_list, method = "t.test") + 
            ggpubr::stat_compare_means(
              method = "anova", 
              label.x = 1,
              label.y = max(.d$p62) + 5, 
              aes(label = paste0(.z, ", Anova, p = ", ..p.format..))
              ) +
            scale_color_brewer(palette = "Set1") +
            theme(
              text = element_text(family = "serif"),
              
              panel.background = element_rect(fill = "white", colour = "black", size = 1),
              axis.title = element_blank(),
              axis.line = element_blank(),
              axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
              
              legend.position = "none"
            ) ->
            .plot
          
          tibble::tibble(
            pval = .pval,
            plot = list(.plot)
          )
        }
        # ggpubr plot
      }
    )
  ) %>% 
  dplyr::select(1, 4) %>% 
  tidyr::unnest() -> 
  p62_subtype

gridExtra::arrangeGrob(
  grobs = p62_subtype %>% dplyr::filter(pval < 0.05) %>% .$plot, 
  nrow = 2, bottom = "Subtype", left = "p62 protein expression"
  ) %>% 
  ggsave(
    filename = "09-p62-cancer-type-subtype.pdf",
    plot = .,
    device = "pdf",
    width = 11,
    height = 10,
    path = path_out
  )


#
#
# end  --------------------------------------------------------------------