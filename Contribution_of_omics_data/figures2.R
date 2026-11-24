# Рисунки 3–6 для отчёта: 3 — итоговый алгоритм, 4 — ассоциации признаков, 5 — генотип × микробиом, 6 — мощность.
# Вход: results/associations.csv, genotype_omics.csv, genotype_omics_input.csv, genotype_permanova.csv,
#       algorithm_models.csv, algorithm_coefs.csv, algorithm_calibration.csv, power.csv
# Запуск из папки analysis:  LC_ALL=en_US.UTF-8 Rscript figures2.R

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(forcats); library(scales); library(patchwork); library(stringr)
})

task_lab <- c(
  CD_vs_UC = "БК против ЯК", IBD_vs_Control = "ВЗК против контроля", Exacerbation = "Клиническое обострение",
  Endoscopic_A2plus = "Эндоскопическая активность A2–A3", Ileal_atrophy = "Атрофия подвздошной кишки",
  CD_ileal = "БК: поражение подвздошной кишки", UC_extensive = "ЯК: тотальный (E3)",
  Advanced_therapy = "Биологическая / JAK-терапия"
)
block_lab <- c("Серология" = "Серология", "Гистология" = "Гистология", "КЦЖК_кал" = "КЦЖК кала",
               "16S" = "Микробиом 16S", "Метаболом_кал" = "Метаболом кала",
               "Метаболом_сыворотка" = "Метаболом сыворотки", "Липидом" = "Липидом",
               "Генетика_WES" = "Генетика", "HLA" = "HLA", "NOD2" = "NOD2")
block_col <- c("Микробиом 16S" = "#1B9E77", "КЦЖК кала" = "#D95F02", "Серология" = "#7570B3",
               "Липидом" = "#E7298A", "Метаболом сыворотки" = "#66A61E", "Метаболом кала" = "#E6AB02",
               "Генетика" = "#A6761D", "HLA" = "#666666", "Гистология" = "#1F78B4", "NOD2" = "#000000")
dec <- function(x, d = 2) formatC(x, format = "f", digits = d, decimal.mark = ",")

theme_report <- theme_minimal(base_size = 10.5) +
  theme(strip.text = element_text(face = "bold", hjust = 0), panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold", size = 12.5), plot.subtitle = element_text(colour = "grey30", size = 9.5),
        plot.caption = element_text(colour = "grey40", size = 8, hjust = 0), legend.position = "top",
        plot.title.position = "plot", plot.caption.position = "plot")
save <- function(p, name, w, h) {
  ggsave(sprintf("results/%s.png", name), p, width = w, height = h, dpi = 300, device = ragg::agg_png, bg = "white")
  ggsave(sprintf("results/%s.pdf", name), p, width = w, height = h, device = cairo_pdf)
}

# короткие имена признаков
short <- function(f) {
  f <- str_replace(f, "^Суммарное содержание С2\\+.*", "Сумма КЦЖК (C2–C6)")
  f <- str_replace(f, "^abs Суммарное содержание ИзоCn.*", "Сумма изо-КЦЖК (абс.)")
  f <- str_replace(f, "^abs (С\\d) \\((.*)\\)", "\\2 кислота (абс.)")
  f <- str_replace(f, "^rel (С\\d) \\((.*)\\)", "\\2 кислота (доля)")
  f <- str_replace(f, "Желчные кислоты в стуле мкмоль/г", "Желчные кислоты кала")
  f <- str_replace(f, "_log2titer", "")
  f <- str_replace_all(f, "_", " ")
  f <- str_replace(f, " \\(1MEOX\\) \\(1TMS\\)", "")
  f <- str_replace(f, "^Pyruvic acid", "Пировиноградная кислота")
  f <- str_replace(f, "NOD2 risk alleles", "NOD2, число аллелей риска")
  f
}

# ---- Рисунок 3: признаки, связанные с клиническими задачами --------------------------------
a <- read.csv("results/associations.csv", check.names = FALSE) %>%
  filter(q < 0.1) %>%
  mutate(task = factor(task_lab[endpoint], levels = task_lab), blk = block_lab[block],
         feat = short(feature),
         robust = case_when(endpoint == "IBD_vs_Control" ~ "Поправка на терапию неприменима (у контролей её нет)",
                            !is.na(q_adj_therapy) & q_adj_therapy < 0.1 ~ "Сохраняется после поправки на терапию",
                            TRUE ~ "Не сохраняется после поправки на терапию")) %>%
  group_by(task) %>% mutate(feat = fct_reorder(paste0(feat, "​", strrep("​", as.integer(task))), effect_sd)) %>% ungroup()

f3 <- ggplot(a, aes(x = effect_sd, y = feat)) +
  geom_vline(xintercept = 0, colour = "grey50") +
  geom_segment(aes(x = 0, xend = effect_sd, yend = feat, colour = blk), linewidth = 0.7) +
  geom_point(aes(colour = blk, shape = robust), size = 2.6, stroke = 1) +
  facet_grid(task ~ ., scales = "free_y", space = "free_y", labeller = label_wrap_gen(18)) +
  scale_colour_manual(values = block_col, name = NULL) +
  scale_shape_manual(values = c("Сохраняется после поправки на терапию" = 16, "Не сохраняется после поправки на терапию" = 1,
                                "Поправка на терапию неприменима (у контролей её нет)" = 17), name = NULL) +
  scale_x_continuous(labels = label_number(accuracy = 0.1, decimal.mark = ",", style_positive = "plus")) +
  labs(title = "Какие показатели отличаются между группами пациентов",
       subtitle = "Признаки с q < 0,1 (поправка Бенджамини–Хохберга внутри блока). Вправо — выше в первой группе задачи (например, при БК, при ВЗК, при обострении)",
       x = "Разница между группами, в стандартных отклонениях (с поправкой на возраст, пол и тип ВЗК)", y = NULL,
       caption = "Визит 1. 16S, метаболомы, липидом — CLR-преобразование долей; КЦЖК — log. Пустой кружок: связь пропадает после поправки на получаемые препараты.") +
  theme_report + theme(strip.text.y = element_text(angle = 0, hjust = 0), legend.box = "vertical",
                       panel.grid.major.y = element_blank())
save(f3, "fig4_feature_associations", 10, 2.2 + 0.2 * nrow(a))

# ---- Рисунок 4: генотип × микробиом ------------------------------------------------------
g <- read.csv("results/genotype_omics_input.csv", check.names = FALSE)
go <- read.csv("results/genotype_omics.csv", check.names = FALSE)
pm <- read.csv("results/genotype_permanova.csv")
panels <- list(
  list(gt = "FUT2_nonsecretor", tx = "Bifidobacterium", lab = c("0" = "секретор", "1" = "non-secretor"), t = "FUT2 → Bifidobacterium"),
  list(gt = "FUT2_nonsecretor", tx = "Faecalibacterium", lab = c("0" = "секретор", "1" = "non-secretor"), t = "FUT2 → Faecalibacterium"),
  list(gt = "NOD2_carrier", tx = "Dialister", lab = c("0" = "нет вариантов", "1" = "носитель"), t = "NOD2 → Dialister"),
  list(gt = "ATG16L1_T300A", tx = "Negativibacillus", lab = c("0" = "0", "1" = "1", "2" = "2"), t = "ATG16L1 T300A (число аллелей) → Negativibacillus")
)
mk <- function(P) {
  d <- g %>% filter(!is.na(.data[[P$gt]])) %>% mutate(G = factor(P$lab[as.character(.data[[P$gt]])], levels = P$lab))
  s <- go %>% filter(genotype == P$gt, feature == P$tx, block == "16S: роды")
  ggplot(d, aes(G, .data[[P$tx]])) +
    geom_boxplot(outlier.shape = NA, fill = "grey92", colour = "grey40", width = 0.55) +
    geom_jitter(aes(colour = diagnosis), width = 0.15, size = 1.4, alpha = 0.8) +
    scale_colour_manual(values = c(CD = "#D95F02", UC = "#1B9E77", Control = "#7570B3"),
                        labels = c(CD = "БК", UC = "ЯК", Control = "Контроль"), name = NULL) +
    labs(title = P$t, subtitle = sprintf("n = %d; эффект %s SD; p = %s; q = %s", nrow(d), dec(s$effect_sd), dec(s$p, 3), dec(s$q)),
         x = NULL, y = "Доля рода (CLR)") +
    theme_report + theme(plot.title = element_text(size = 10.5), plot.subtitle = element_text(size = 8.5))
}
pm_txt <- paste0("Общий состав микробиома (PERMANOVA, с поправкой на диагноз, возраст, пол): ",
                 paste(sprintf("%s — R² = %s%%, p = %s", c("FUT2", "NOD2", "ATG16L1"), dec(100 * pm$R2, 1), dec(pm$p, 2)), collapse = "; "), ".")
f4 <- wrap_plots(lapply(panels, mk), ncol = 2, guides = "collect") +
  plot_annotation(title = "Генотип и микробиом кишечника",
                  subtitle = "Роды бактерий с наибольшей разницей между генотипами (16S, визит 1); модель с поправкой на диагноз, возраст и пол",
                  caption = paste0(pm_txt, "\nnon-secretor FUT2 — только 5 человек, поэтому эффекты для FUT2 оценены неточно."),
                  theme = theme_report) & theme(legend.position = "bottom")
save(f4, "fig5_genotype_microbiome", 10, 8)

# ---- Рисунок 5: итоговый алгоритм ---------------------------------------------------------
m <- read.csv("results/algorithm_models.csv") %>% mutate(task = factor(task_lab[endpoint], levels = rev(task_lab)))
mm <- bind_rows(
  m %>% transmute(task, model = "Базовая модель (клиника + рутинные анализы)", auc, lo = auc_lo, hi = auc_hi),
  m %>% filter(auc_calpro >= 0.5) %>%   # AUC < 0,5 вне выборки = сигнала нет, оценка «перевёрнута»
    transmute(task, model = "Только кальпротектин", auc = auc_calpro, lo = auc_calpro_lo, hi = auc_calpro_hi))
f5a <- ggplot(mm, aes(auc, task, colour = model)) +
  annotate("rect", xmin = -Inf, xmax = 0.5, ymin = -Inf, ymax = Inf, fill = "grey93") +
  geom_vline(xintercept = c(0.5, 0.7, 0.8), linetype = c("solid", "dotted", "dotted"), colour = "grey55") +
  geom_errorbar(aes(xmin = lo, xmax = hi), width = 0, linewidth = 0.8, orientation = "y", position = position_dodge(0.55)) +
  geom_point(size = 2.5, position = position_dodge(0.55)) +
  scale_colour_manual(values = c("Базовая модель (клиника + рутинные анализы)" = "#2166AC", "Только кальпротектин" = "#B2182B"), name = NULL) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1), labels = label_number(accuracy = 0.1, decimal.mark = ",")) +
  labs(title = "А. Качество предсказания на новых пациентах", x = "AUC (95% ДИ); серая зона — хуже угадывания", y = NULL,
       caption = "Кальпротектин не показан там, где сам по себе он не предсказывает (AUC < 0,5).") +
  theme_report + theme(panel.grid.major.y = element_blank(), legend.direction = "vertical")

cf <- read.csv("results/algorithm_coefs.csv") %>%
  filter(endpoint %in% c("IBD_vs_Control", "Endoscopic_A2plus", "Exacerbation", "Advanced_therapy")) %>%
  mutate(task = factor(task_lab[endpoint], levels = task_lab), sig = p < 0.05,
         label = factor(label, levels = rev(unique(label))))
f5b <- ggplot(cf, aes(OR, label, colour = sig)) +
  geom_vline(xintercept = 1, colour = "grey50") +
  geom_errorbar(aes(xmin = OR_lo, xmax = OR_hi), width = 0, orientation = "y") + geom_point(size = 2) +
  facet_wrap(~task, ncol = 2, scales = "free_y") +
  scale_x_log10(labels = label_number(accuracy = 0.1, decimal.mark = ",")) +
  coord_cartesian(xlim = c(0.1, 30)) +
  scale_colour_manual(values = c(`TRUE` = "#2166AC", `FALSE` = "grey60"), labels = c(`TRUE` = "p < 0,05", `FALSE` = "p ≥ 0,05"), name = NULL) +
  labs(title = "В. Вклад показателей в базовую модель", x = "Отношение шансов на 1 стандартное отклонение (95% ДИ, лог. шкала)", y = NULL) +
  theme_report

cal <- read.csv("results/algorithm_calibration.csv") %>%
  filter(endpoint %in% c("IBD_vs_Control", "Endoscopic_A2plus", "Advanced_therapy")) %>%
  mutate(task = factor(task_lab[endpoint], levels = task_lab))
f5c <- ggplot(cal, aes(pred, obs, colour = task)) +
  geom_abline(linetype = "dashed", colour = "grey50") + geom_line() + geom_point(aes(size = n)) +
  scale_size_continuous(range = c(1.2, 3), guide = "none") +
  scale_colour_brewer(palette = "Dark2", name = NULL) +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  scale_x_continuous(labels = label_percent()) + scale_y_continuous(labels = label_percent()) +
  labs(title = "Б. Калибровка: совпадает ли предсказанный риск с наблюдаемым", x = "Предсказанный риск (квинтили)", y = "Наблюдаемая частота") +
  theme_report + theme(legend.direction = "vertical")

f5 <- (f5a | f5c) / f5b + plot_layout(heights = c(1, 1.15)) +
  plot_annotation(title = "Итоговый алгоритм на рутинных данных",
                  caption = "Визит 1. AUC — 20 повторов 5-кратной кросс-валидации, ДИ — бутстреп. Коэффициенты — логистическая регрессия на всех пациентах задачи, пропуски заменены медианой.",
                  theme = theme_report)
save(f5, "fig3_algorithm", 11, 10)

# ---- Рисунок 6: мощность --------------------------------------------------------------------
pw <- read.csv("results/power.csv") %>% mutate(blk = block_lab[block])
f6 <- ggplot(pw, aes(n, mde, colour = blk)) +
  geom_hline(yintercept = 0.05, linetype = "dashed", colour = "grey40") +
  annotate("text", x = max(pw$n), y = 0.053, label = "прирост +0,05 AUC", hjust = 1, vjust = 0, size = 3, colour = "grey30") +
  geom_point(size = 2.3, alpha = 0.85) +
  scale_colour_manual(values = block_col, name = NULL) +
  scale_y_continuous(labels = label_number(accuracy = 0.01, decimal.mark = ",", style_positive = "plus"), limits = c(0, NA)) +
  labs(title = "Какой прирост можно было обнаружить при таком числе пациентов",
       subtitle = "Минимальный прирост AUC, который выявляется с вероятностью 80% (α = 0,05); каждая точка — пара «задача × блок»",
       x = "Число пациентов с данным блоком", y = "Минимальный обнаруживаемый ΔAUC",
       caption = "Оценено по перестановочному распределению прироста (случайный блок на тех же пациентах): 95-й перцентиль + 0,84 × SD.\nТочки выше пунктира: реальный прирост +0,05 в этой паре, скорее всего, остался бы незамеченным.") +
  theme_report + theme(legend.position = "right")
save(f6, "fig6_power", 9, 5.5)

cat("Готово: fig3–fig6 (+ PDF)\n")
