# Рисунки для отчёта: прирост AUC от омиксных блоков к базовой модели.
# Вход:  results/increment_results.json (создаётся run_increment.py)
# Выход: results/fig1_delta_auc.{png,pdf}, results/fig2_auc_base_vs_extended.{png,pdf}
# Запуск из папки analysis:  LC_ALL=en_US.UTF-8 Rscript figures.R
# (в RStudio — просто Source; файл в кодировке UTF-8)

suppressPackageStartupMessages({
  library(jsonlite); library(dplyr); library(ggplot2); library(forcats); library(scales)
})

res <- fromJSON("results/increment_results.json")
d <- as_tibble(res$results)

# --- подписи по-русски -------------------------------------------------------
task_lab <- c(
  CD_vs_UC          = "Болезнь Крона против язвенного колита",
  IBD_vs_Control    = "ВЗК против контроля",
  Exacerbation      = "Клиническое обострение",
  Endoscopic_A2plus = "Эндоскопическая активность A2–A3",
  Ileal_atrophy     = "Атрофия подвздошной кишки"
)
block_lab <- c(
  "Серология" = "Серология", "Гистология" = "Гистология",
  "КЦЖК_кал" = "КЦЖК кала", "16S" = "Микробиом 16S",
  "Метаболом_кал" = "Метаболом кала", "Метаболом_сыворотка" = "Метаболом сыворотки",
  "Липидом" = "Липидом"
)

d <- d %>%
  mutate(
    ci_lo   = sapply(delta_ci, `[`, 1), ci_hi = sapply(delta_ci, `[`, 2),
    q       = p.adjust(p_perm, method = "BH"),
    task    = factor(task_lab[endpoint], levels = task_lab),
    block_r = block_lab[block],
    label   = sprintf("%s (n = %d)", block_r, n),
    # кандидат для проверки: ΔAUC ≥ 0,02 и p < 0,15, и блок сам по себе предсказывает
    # в нужную сторону (AUC блока ≥ 0,5; иначе прирост — артефакт стекинга)
    status  = case_when(delta >= 0.02 & p_perm < 0.15 & auc_block >= 0.5 ~ "Кандидат для проверки (ΔAUC ≥ 0,02, p < 0,15)",
                        TRUE ~ "Прироста нет"),
    p_txt   = ifelse(p_perm < 0.01, "p < 0,01", paste0("p = ", formatC(p_perm, format = "f", digits = 2, decimal.mark = ",")))
  ) %>%
  group_by(task) %>% mutate(label = fct_reorder(label, delta)) %>% ungroup()

col_main <- c("Кандидат для проверки (ΔAUC ≥ 0,02, p < 0,15)" = "#2166AC", "Прироста нет" = "grey55")

theme_report <- theme_minimal(base_size = 11) +
  theme(
    strip.text = element_text(face = "bold", hjust = 0, size = 10.5),
    panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(colour = "grey30", size = 10),
    plot.caption = element_text(colour = "grey40", size = 8.5, hjust = 0),
    legend.position = "top", legend.title = element_blank(),
    plot.title.position = "plot", plot.caption.position = "plot"
  )

# --- Рисунок 1: прирост AUC с 95% ДИ ----------------------------------------
x_max <- max(d$ci_hi) + 0.06
f1 <- ggplot(d, aes(y = label, x = delta, colour = status)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40") +
  geom_errorbar(aes(xmin = ci_lo, xmax = ci_hi), width = 0, linewidth = 0.8, orientation = "y") +
  geom_point(size = 2.4) +
  geom_text(aes(x = x_max, label = p_txt), hjust = 1, size = 3, colour = "grey35") +
  facet_grid(task ~ ., scales = "free_y", space = "free_y", switch = "y",
             labeller = label_wrap_gen(width = 22)) +
  scale_colour_manual(values = col_main) +
  scale_x_continuous(labels = label_number(accuracy = 0.01, decimal.mark = ",", style_positive = "plus"),
                     breaks = seq(-0.1, 0.15, 0.05)) +
  labs(
    title = "Омиксные данные почти не улучшают предсказание по сравнению с рутинными анализами",
    subtitle = "Прирост AUC (ΔAUC) при добавлении блока к базовой модели (клиника + рутинные анализы); точка — среднее, линия — 95% ДИ",
    x = "ΔAUC = AUC(база + блок) − AUC(база)", y = NULL,
    caption = "Визит 1. 20 повторов 5-кратной кросс-валидации; ДИ — бутстреп 2000 выборок; p — перестановочный тест (100 перестановок).\nПосле поправки Бенджамини–Хохберга на 29 сравнений ни один прирост не достоверен (минимальное q = 0,29)."
  ) +
  theme_report +
  theme(strip.placement = "outside", strip.text.y.left = element_text(angle = 0, vjust = 1))

ggsave("results/fig1_delta_auc.png", f1, width = 10, height = 9.5, dpi = 300, device = ragg::agg_png, bg = "white")
ggsave("results/fig1_delta_auc.pdf", f1, width = 10, height = 9.5, device = cairo_pdf)

# --- Рисунок 2: AUC базовой и расширенной модели ----------------------------
d2 <- d %>% mutate(label = fct_reorder(label, auc_m0))
f2 <- ggplot(d2, aes(y = label)) +
  geom_vline(xintercept = 0.5, linetype = "dotted", colour = "grey50") +
  geom_segment(aes(x = auc_m0, xend = auc_m1, yend = label), colour = "grey70", linewidth = 1) +
  geom_point(aes(x = auc_m0, shape = "База: клиника + рутинные анализы"), size = 2.6, colour = "grey30") +
  geom_point(aes(x = auc_m1, shape = "База + блок"), size = 2.6, colour = "#2166AC") +
  geom_text(data = tibble(task = factor(task_lab[1], levels = task_lab)), aes(x = 0.5, y = Inf, label = "угадывание"),
            inherit.aes = FALSE, vjust = 1.5, hjust = -0.1, size = 3, colour = "grey45") +
  facet_grid(task ~ ., scales = "free_y", space = "free_y", switch = "y",
             labeller = label_wrap_gen(width = 22)) +
  scale_shape_manual(values = c("База: клиника + рутинные анализы" = 16, "База + блок" = 17), breaks = c("База: клиника + рутинные анализы", "База + блок")) +
  scale_x_continuous(limits = c(0.4, 1), breaks = seq(0.4, 1, 0.1),
                     labels = label_number(accuracy = 0.1, decimal.mark = ",")) +
  labs(
    title = "Базовая модель уже достигает почти того же качества, что и расширенная",
    subtitle = "AUC базовой модели (серый круг) и базы с добавленным блоком (синий треугольник) на одних и тех же пациентах",
    x = "AUC (0,5 — угадывание, 1,0 — идеальное разделение)", y = NULL,
    caption = "Визит 1. Оценки получены на пациентах, не участвовавших в обучении (повторная кросс-валидация)."
  ) +
  theme_report +
  theme(strip.placement = "outside", strip.text.y.left = element_text(angle = 0, vjust = 1))

ggsave("results/fig2_auc_base_vs_extended.png", f2, width = 10, height = 9.5, dpi = 300, device = ragg::agg_png, bg = "white")
ggsave("results/fig2_auc_base_vs_extended.pdf", f2, width = 10, height = 9.5, device = cairo_pdf)

cat("Готово: results/fig1_delta_auc.png, results/fig2_auc_base_vs_extended.png (+ PDF)\n")
