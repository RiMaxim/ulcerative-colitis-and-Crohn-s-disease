# ==============================================================================
# ОПТИМИЗИРОВАННЫЙ КОД
# ==============================================================================

library(tidyverse)
library(tibble)
library(pROC)  # для расчета AUC
library(purrr)  # для map функций

# ==============================================================================
# ЗАГРУЗКА ДАННЫХ
# ==============================================================================

emu <- read.table("Desktop/WORK/gut/input_SCFA.txt", 
                  header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)

emu <- emu %>% rownames_to_column("sample")

emu_long <- emu %>% 
  pivot_longer(cols = -sample, names_to = "species", values_to = "abundance") %>%
  mutate(
    species = gsub("_", " ", species),
    species = gsub("\\[|\\]", "", species),
    species = gsub("-", " ", species),
    species = tolower(trimws(species))
  )

# ==============================================================================
# ЗАГРУЗКА РЕФЕРЕНСНЫХ БАЗ
# ==============================================================================

ref_db1 <- read.table("Desktop/WORK/gut/Propionate.txt", sep = "\t", header = TRUE, check.names = FALSE)
ref_db2 <- read.table("Desktop/WORK/gut/Butyrate.txt", sep = "\t", header = TRUE, check.names = FALSE)
ref_db3 <- read.table("Desktop/WORK/gut/AFL.txt", sep = "\t", header = TRUE, check.names = FALSE)

# ==============================================================================
# ФУНКЦИЯ ДЛЯ ОБРАБОТКИ РЕФЕРЕНСНЫХ ДАННЫХ
# ==============================================================================

process_reference <- function(data, phenotype_col) {
  data %>%
    filter(!is.na(Species)) %>%
    filter(.data[[phenotype_col]] == 1) %>%
    mutate(
      species = gsub("_", " ", Species),
      species = gsub("-", " ", species),
      species = gsub("\\[|\\]", "", species),
      species = tolower(trimws(species))
    ) %>%
    select(species) %>%
    distinct()
}

# Обрабатываем все референсы
ref_propionate <- process_reference(ref_db1, "Binary Phenotype")
ref_butyrate <- process_reference(ref_db2, "Binary Phenotype")
ref_acetate <- process_reference(ref_db3, "Acetate Binary Phenotype")
ref_formate <- process_reference(ref_db3, "Formate Binary Phenotype")
ref_l_lactate <- process_reference(ref_db3, "L-Lactate Binary Phenotype")
ref_d_lactate <- process_reference(ref_db3, "D-Lactate Binary Phenotype")

# ==============================================================================
# ФУНКЦИЯ ДЛЯ РАСЧЁТА CPI
# ==============================================================================

calculate_cpi <- function(data, ref_data, index_name) {
  data %>%
    filter(species %in% ref_data$species) %>%
    group_by(sample) %>%
    summarise(
      !!index_name := sum(abundance, na.rm = TRUE),
      .groups = "drop"
    )
}

# Рассчитываем все CPI
cpi_propionate <- calculate_cpi(emu_long, ref_propionate, "CPI_Propionate")
cpi_butyrate <- calculate_cpi(emu_long, ref_butyrate, "CPI_Butyrate")
cpi_acetate <- calculate_cpi(emu_long, ref_acetate, "CPI_Acetate")
cpi_formate <- calculate_cpi(emu_long, ref_formate, "CPI_Formate")
cpi_l_lactate <- calculate_cpi(emu_long, ref_l_lactate, "CPI_L_Lactate")
cpi_d_lactate <- calculate_cpi(emu_long, ref_d_lactate, "CPI_D_Lactate")

# ==============================================================================
# ОБЪЕДИНЕНИЕ ВСЕХ CPI В ОДНУ ТАБЛИЦУ
# ==============================================================================

community_phenotype_index_all <- cpi_propionate %>%
  full_join(cpi_butyrate, by = "sample") %>%
  full_join(cpi_acetate, by = "sample") %>%
  full_join(cpi_formate, by = "sample") %>%
  full_join(cpi_l_lactate, by = "sample") %>%
  full_join(cpi_d_lactate, by = "sample") %>%
  mutate(across(starts_with("CPI_"), ~replace_na(.x, 0)))

# ==============================================================================
# СТАТИСТИКА ПОКРЫТИЯ (ВСЕ В ОДНОЙ ТАБЛИЦЕ)
# ==============================================================================

coverage_stats <- data.frame(
  Phenotype = c("Propionate", "Butyrate", "Acetate", "Formate", "L-Lactate", "D-Lactate"),
  Matched = c(
    length(intersect(emu_long$species, ref_propionate$species)),
    length(intersect(emu_long$species, ref_butyrate$species)),
    length(intersect(emu_long$species, ref_acetate$species)),
    length(intersect(emu_long$species, ref_formate$species)),
    length(intersect(emu_long$species, ref_l_lactate$species)),
    length(intersect(emu_long$species, ref_d_lactate$species))
  ),
  Total_EMU = length(unique(emu_long$species))
) %>%
  mutate(
    Percent = round(100 * Matched / Total_EMU, 2)
  )

print("=== SPECIES COVERAGE STATISTICS ===")
print(coverage_stats)

# ==============================================================================
# СОХРАНЕНИЕ РЕЗУЛЬТАТОВ
# ==============================================================================

write.table(coverage_stats, "coverage_statistics.tsv", 
            sep = "\t", row.names = FALSE, quote = FALSE)

write.table(community_phenotype_index_all, "CPI_all_phenotypes.tsv", 
            sep = "\t", row.names = FALSE, quote = FALSE)

# ==============================================================================
# ПРОСМОТР РЕЗУЛЬТАТОВ
# ==============================================================================

print("=== COMMUNITY PHENOTYPE INDICES ===")
print(head(community_phenotype_index_all))














# ==============================================================================
# ПОЛНЫЙ КОД ДЛЯ ВИЗУАЛИЗАЦИИ CPI С P-VALUE
# ==============================================================================

library(tidyverse)
library(ggplot2)
library(ggpubr)
library(patchwork)
library(RColorBrewer)
library(rstatix)

# ==============================================================================
# 1. ИЗВЛЕЧЕНИЕ МЕТАДАННЫХ
# ==============================================================================

cpi_with_metadata <- community_phenotype_index_all %>%
  separate(sample, into = c("ID1", "Patient_ID", "Visit", "Disease", "Sex", "Age", "ID2"), 
           sep = "\\|", remove = FALSE, fill = "right") %>%
  mutate(
    Patient_ID = as.integer(Patient_ID),
    Visit = as.integer(Visit),
    Age = suppressWarnings(as.integer(Age)),
    Disease = case_when(
      Disease == "CD" ~ "Crohn's Disease",
      Disease == "UC" ~ "Ulcerative Colitis",
      Disease == "Control" ~ "Healthy Control",
      is.na(Disease) | Disease == "" ~ "Unknown",
      TRUE ~ Disease
    ),
    Sex = case_when(
      Sex == "M" ~ "Male",
      Sex == "F" ~ "Female",
      is.na(Sex) | Sex == "" ~ "Unknown",
      TRUE ~ Sex
    ),
    Age_Group = case_when(
      Age <= 18 ~ "0-18",
      Age <= 40 ~ "19-40",
      Age <= 60 ~ "41-60",
      Age > 60 ~ "60+",
      TRUE ~ "Unknown"
    ),
    Disease = factor(Disease, levels = c("Healthy Control", "Crohn's Disease", 
                                         "Ulcerative Colitis", "Unknown")),
    Sex = factor(Sex, levels = c("Male", "Female", "Unknown")),
    Age_Group = factor(Age_Group, levels = c("0-18", "19-40", "41-60", "60+", "Unknown"))
  ) %>%
  select(-ID1, -ID2)

# ==============================================================================
# 2. ПЕРЕВОД В LONG FORMAT
# ==============================================================================

cpi_disease_long <- cpi_with_metadata %>%
  pivot_longer(
    cols = starts_with("CPI_"),
    names_to = "Phenotype",
    values_to = "CPI"
  ) %>%
  mutate(
    Phenotype = gsub("CPI_", "", Phenotype),
    Phenotype = factor(Phenotype, levels = c("Propionate", "Butyrate", "Acetate", 
                                             "Formate", "L_Lactate", "D_Lactate"))
  ) %>%
  filter(!is.na(CPI))

# ==============================================================================
# 3 и 4. ИСПРАВЛЕННЫЙ РАСЧЁТ P-VALUE И ТОЧНЫХ КООРДИНАТ ДЛЯ КАЖДОГО ФАСЕТА
# ==============================================================================

# Общая функция для кастомизации подписей, цветов и шрифтов
format_labels <- function(data) {
  data %>%
    mutate(
      p_label = case_when(
        p.adj < 0.001 ~ "Padj < 0.001***",
        p.adj < 0.01 ~ paste0("Padj = ", round(p.adj, 3), "**"),
        p.adj < 0.05 ~ paste0("Padj = ", round(p.adj, 3), "*"),
        TRUE ~ paste0("Padj = ", round(p.adj, 3))
      ),
      is_significant = p.adj < 0.05,
      color = ifelse(is_significant, "#D55E00", "#7f8c8d"), 
      face = ifelse(is_significant, "bold", "plain")
    )
}

# --- ДЛЯ DISEASE ---
global_max_d <- max(p1_data$CPI, na.rm = TRUE)
step_d <- global_max_d * 0.12

disease_p_values <- p1_data %>%
  group_by(Phenotype) %>%
  wilcox_test(CPI ~ Disease) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance() %>%
  mutate(
    x = as.numeric(factor(group1, levels = levels(p1_data$Disease))),
    xend = as.numeric(factor(group2, levels = levels(p1_data$Disease))),
    pair_id = paste0(group1, "_vs_", group2),
    y.position = case_when(
      pair_id == "Healthy Control_vs_Crohn's Disease"   ~ global_max_d + step_d * 0.4,
      pair_id == "Crohn's Disease_vs_Ulcerative Colitis" ~ global_max_d + step_d * 0.4,
      pair_id == "Healthy Control_vs_Ulcerative Colitis" ~ global_max_d + step_d * 1.5,
      TRUE ~ global_max_d + step_d
    )
  ) %>% format_labels()

# --- ДЛЯ SEX ---
global_max_s <- max(p2_data$CPI, na.rm = TRUE)
step_s <- global_max_s * 0.12

sex_p_values <- p2_data %>%
  group_by(Phenotype) %>%
  wilcox_test(CPI ~ Sex) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance() %>%
  mutate(
    x = as.numeric(factor(group1, levels = levels(p2_data$Sex))),
    xend = as.numeric(factor(group2, levels = levels(p2_data$Sex))),
    y.position = global_max_s + step_s * 0.4
  ) %>% format_labels()

# --- ДЛЯ AGE GROUP (ИСПРАВЛЕНО РАНЖИРОВАНИЕ ЯРУСОВ) ---
global_max_a <- max(p3_data$CPI, na.rm = TRUE)
step_a <- global_max_a * 0.10 # Слегка уменьшили шаг для компактности

age_p_values <- p3_data %>%
  group_by(Phenotype) %>%
  wilcox_test(CPI ~ Age_Group) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance() %>%
  mutate(
    x = as.numeric(factor(group1, levels = levels(p3_data$Age_Group))),
    xend = as.numeric(factor(group2, levels = levels(p3_data$Age_Group))),
    dist = abs(xend - x)
  ) %>% 
  # ВАЖНОЕ ИСПРАВЛЕНИЕ: группируем только по Phenotype, 
  # а ярус (tier) определяем по плотному рангу дистанции внутри этого метаболита
  group_by(Phenotype) %>% 
  mutate(tier = dense_rank(dist)) %>% 
  ungroup() %>%
  mutate(
    # Теперь ярусы гарантированно будут идти строго по порядку: 1, 2, 3... без пропусков
    y.position = global_max_a + (step_a * 0.3) + (step_a * 0.8 * (tier - 1))
  ) %>% 
  format_labels()



# ==============================================================================
# 5. ПОДГОТОВКА ДАННЫХ ДЛЯ ГРАФИКОВ
# ==============================================================================

# Данные для Disease
p1_data <- cpi_disease_long %>%
  filter(Disease != "Unknown")

# Данные для Sex
p2_data <- cpi_disease_long %>%
  filter(Sex != "Unknown")

# Данные для Age Group
p3_data <- cpi_disease_long %>%
  filter(Age_Group != "Unknown")

# ==============================================================================
# 6. СБОРКА ГРАФИКА (С минимальными отступами осей)
# ==============================================================================

# --- ГРАФИК 1 (DISEASE) ---
p1 <- ggplot(p1_data, aes(x = Disease, y = CPI, fill = Disease)) +
  geom_boxplot(alpha = 0.7, outlier.size = 1.5) +
  geom_jitter(width = 0.2, alpha = 0.4, size = 1) +
  facet_wrap(~Phenotype, scales = "fixed", ncol = 3) + # ИСПРАВЛЕНО: ось зафиксирована
  scale_fill_brewer(palette = "Set2") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 12),
    legend.position = "none",
    strip.text = element_text(size = 11, face = "bold"),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    plot.margin = margin(t = 20, r = 10, b = 10, l = 10)
  ) +
  labs(title = "Community Phenotype Index by Disease Type", x = "", y = "CPI value")

p1_with_p <- p1 +
  stat_pvalue_manual(disease_p_values, label = "p_label", vjust = -0.5, size = 2.8,
                     tip.length = 0.02, bracket.shorten = 0.05, step.increase = 0,
                     hide.ns = FALSE, color = "color", fontface = "face") + 
  scale_y_continuous(
    limits = c(0, 1.3),                    # Жесткие границы от 0 до 4.5
    breaks = seq(0, 1.3, by = 0.5),        # Шаг сетки на оси через каждые 0.5
    expand = expansion(mult = c(0.05, 0))   # 0% отступа снизу (ровно от нуля) и 10% запаса сверху
  )

# --- ГРАФИК 2 (SEX) ---
p2 <- ggplot(p2_data, aes(x = Sex, y = CPI, fill = Sex)) +
  geom_boxplot(alpha = 0.7, outlier.size = 1.5) +
  geom_jitter(width = 0.2, alpha = 0.4, size = 1) +
  facet_wrap(~Phenotype, scales = "fixed", ncol = 3) + # ИСПРАВЛЕНО: ось зафиксирована
  scale_fill_manual(values = c("Male" = "#4E79A7", "Female" = "#F28E2B")) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 12),
    legend.position = "none",
    strip.text = element_text(size = 11, face = "bold"),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    plot.margin = margin(t = 20, r = 10, b = 10, l = 10)
  ) +
  labs(title = "Community Phenotype Index by Sex", x = "", y = "CPI value")

p2_with_p <- p2 +
  stat_pvalue_manual(sex_p_values, label = "p_label", vjust = -0.5, size = 2.8,
                     tip.length = 0.02, bracket.shorten = 0.05, step.increase = 0,
                     hide.ns = FALSE, color = "color", fontface = "face") + 
  scale_y_continuous(
    limits = c(0, 1.15),                    # Жесткие границы от 0 до 4.5
    breaks = seq(0, 1.15, by = 0.5),        # Шаг сетки на оси через каждые 0.5
    expand = expansion(mult = c(0.05, 0))   # 0% отступа снизу (ровно от нуля) и 10% запаса сверху
  )


# --- ГРАФИК 3 (AGE GROUP) ---
p3 <- ggplot(p3_data, aes(x = Age_Group, y = CPI, fill = Age_Group)) +
  geom_boxplot(alpha = 0.7, outlier.size = 1.5) +
  geom_jitter(width = 0.2, alpha = 0.4, size = 1) +
  facet_wrap(~Phenotype, scales = "fixed", ncol = 3) + # ИСПРАВЛЕНО: ось зафиксирована
  scale_fill_brewer(palette = "RdYlBu") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 12),
    legend.position = "none",
    strip.text = element_text(size = 11, face = "bold"),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    plot.margin = margin(t = 20, r = 10, b = 10, l = 10)
  ) +
  labs(title = "Community Phenotype Index by Age Group", x = "", y = "CPI value")

p3_with_p <- p3 +
  stat_pvalue_manual(age_p_values, label = "p_label", vjust = -0.5, size = 2.8,
                     tip.length = 0.02, bracket.shorten = 0.05, step.increase = 0,
                     hide.ns = FALSE, color = "color", fontface = "face") +
  scale_y_continuous(
    limits = c(0, 1.2),                    # Жесткие границы от 0 до 4.5
    breaks = seq(0, 1.2, by = 0.5),        # Шаг сетки на оси через каждые 0.5
    expand = expansion(mult = c(0.05, 0))   # 0% отступа снизу (ровно от нуля) и 10% запаса сверху
  )

# Сохранение результатов
ggsave("Desktop/1.pdf", p1_with_p, width = 14, height = 8)
ggsave("Desktop/2.pdf", p2_with_p, width = 14, height = 8)
ggsave("Desktop/3.pdf", p3_with_p, width = 14, height = 8)


# ==============================================================================
# 9. ГРАФИК 4: ТЕПЛОВАЯ КАРТА
# ==============================================================================

mean_cpi <- cpi_disease_long %>%
  filter(Disease != "Unknown", Sex != "Unknown") %>%
  group_by(Disease, Sex, Phenotype) %>%
  summarise(
    mean_CPI = mean(CPI, na.rm = TRUE),
    count = n(),
    .groups = "drop"
  )

p4 <- ggplot(mean_cpi, aes(x = Disease, y = Phenotype, fill = mean_CPI)) +
  geom_tile(color = "white", linewidth = 0.5) +
  facet_wrap(~Sex) +
  scale_fill_gradient2(
    low = "#2166AC", 
    mid = "white", 
    high = "#B2182B",
    midpoint = median(mean_cpi$mean_CPI, na.rm = TRUE),
    name = "Mean CPI"
  ) +
  geom_text(aes(label = round(mean_CPI, 3)), size = 3.5, fontface = "bold") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 12),
    axis.text.y = element_text(size = 12),
    strip.text = element_text(size = 11, face = "bold"),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  ) +
  labs(
    title = "Heatmap of Mean CPI by Disease and Sex",
    x = "",
    y = ""
  )

ggsave("Desktop/4.pdf", p4, width = 14, height = 8)



# ==============================================================================
# ГРАФИК: ВЛИЯНИЕ ВОЗРАСТА НА ВСЕ ФЕНОТИПЫ CPI С РАЗДЕЛЕНИЕМ ПО БОЛЕЗНЯМ
# ==============================================================================

p5 <- analysis_data %>%
  ggplot(aes(x = Age, y = CPI, color = Disease)) +  # ← изменили color с Phenotype на Disease
  geom_point(alpha = 0.2, size = 1) +
  geom_smooth(method = "loess", se = TRUE, span = 0.8) +  # теперь кривые для каждой болезни
  facet_wrap(~Phenotype, scales = "free_y") +  # отдельные панели для каждого фенотипа
  scale_color_brewer(palette = "Set1", name = "Disease") +  # яркие цвета для болезней
  theme_minimal() +
  labs(title = "Age Effects on All CPI Phenotypes by Disease",
       subtitle = "LOESS smoothing with 95% CI",
       x = "Age (years)", y = "CPI Value") +
  theme(
    legend.position = "bottom",  # показываем легенду внизу
    strip.text = element_text(face = "bold", size = 11),
    legend.title = element_text(face = "bold"),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12, face = "bold")
  )

# Сохранение
ggsave("Desktop/5.pdf", p5, width = 14, height = 8)



# ==============================================================================
# АНАЛИЗ ДИНАМИКИ CPI МЕЖДУ ВИЗИТАМИ
# ==============================================================================

library(tidyverse)
library(ggplot2)
library(patchwork)

# ==============================================================================
# 1. ПОДГОТОВКА ДАННЫХ ДЛЯ АНАЛИЗА ВИЗИТОВ
# ==============================================================================

# Преобразуем данные в формат с визитами
visit_data <- cpi_with_metadata %>%
  # Оставляем только пациентов с несколькими визитами
  group_by(Patient_ID) %>%
  filter(n() > 1) %>%  # только те, у кого >1 визита
  ungroup() %>%
  # Переводим в long format
  pivot_longer(
    cols = starts_with("CPI_"),
    names_to = "Phenotype",
    values_to = "CPI"
  ) %>%
  mutate(
    Phenotype = gsub("CPI_", "", Phenotype),
    Visit = factor(Visit, levels = c("1", "2", "3", "4"))  # упорядочиваем визиты
  ) %>%
  filter(!is.na(CPI), Disease != "Unknown")

# ==============================================================================
# 2. ГРАФИК ИНДИВИДУАЛЬНЫХ ТРАЕКТОРИЙ
# ==============================================================================

# Для каждого пациента показываем, как меняется CPI между визитами
p_individual <- visit_data %>%
  filter(Phenotype == "Propionate") %>%  # выберите конкретный фенотип
  ggplot(aes(x = Visit, y = CPI, group = Patient_ID, color = Disease)) +
  geom_line(alpha = 0.5, size = 0.8) +
  geom_point(alpha = 0.6, size = 1.5) +
  facet_wrap(~Disease, ncol = 3) +
  scale_color_brewer(palette = "Set1") +
  theme_minimal() +
  labs(
    title = "Individual Trajectories of Propionate CPI Across Visits",
    subtitle = "Each line represents one patient",
    x = "Visit Number",
    y = "CPI Propionate"
  ) +
  theme(
    legend.position = "bottom",
    strip.text = element_text(face = "bold", size = 11),
    axis.text = element_text(size = 10)
  )

ggsave("Desktop/individual_trajectories.pdf", p_individual, width = 12, height = 8)

# ==============================================================================
# 3. СРЕДНИЕ ТРАЕКТОРИИ ПО ГРУППАМ
# ==============================================================================

# Средние значения и стандартные ошибки для каждой группы
mean_trajectories <- visit_data %>%
  group_by(Disease, Phenotype, Visit) %>%
  summarise(
    mean_CPI = mean(CPI, na.rm = TRUE),
    se_CPI = sd(CPI, na.rm = TRUE) / sqrt(n()),
    n = n(),
    .groups = "drop"
  )

# График средних траекторий
p_mean <- mean_trajectories %>%
  filter(Phenotype == "Propionate") %>%
  ggplot(aes(x = Visit, y = mean_CPI, color = Disease, group = Disease)) +
  geom_line(size = 1.2) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = mean_CPI - se_CPI, ymax = mean_CPI + se_CPI),
                width = 0.2, alpha = 0.5) +
  scale_color_brewer(palette = "Set1", name = "Disease") +
  theme_minimal() +
  labs(
    title = "Mean Propionate CPI Trajectories by Disease",
    subtitle = "Points = mean ± SE",
    x = "Visit Number",
    y = "Mean CPI Propionate"
  ) +
  theme(
    legend.position = "bottom",
    axis.text = element_text(size = 11),
    axis.title = element_text(size = 12, face = "bold")
  )

ggsave("Desktop/mean_trajectories.pdf", p_mean, width = 10, height = 6)

# ==============================================================================
# 4. ВСЕ ФЕНОТИПЫ В ОДНОЙ ПАНЕЛИ
# ==============================================================================

p6 <- mean_trajectories %>%
  ggplot(aes(x = Visit, y = mean_CPI, color = Disease, group = Disease)) +
  geom_line(size = 1.2) +
  geom_point(size = 2.5) +
  geom_errorbar(aes(ymin = mean_CPI - se_CPI, ymax = mean_CPI + se_CPI),
                width = 0.2, alpha = 0.4) +
  facet_wrap(~Phenotype, scales = "free_y", ncol = 3) +
  scale_color_brewer(palette = "Set1", name = "Disease") +
  theme_minimal() +
  labs(
    title = "Mean CPI Trajectories by Phenotype and Disease",
    subtitle = "Points = mean ± SE",
    x = "Visit Number",
    y = "Mean CPI"
  ) +
  theme(
    legend.position = "bottom",
    strip.text = element_text(face = "bold", size = 11),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12, face = "bold")
  )

ggsave("Desktop/6.pdf", p6, width = 14, height = 8)


#############################Combine

combined_plot <- (p1_with_p + p2_with_p) / (p3_with_p + p5) / (p4 + p6) +
  plot_annotation (title = "Community Phenotype Index Analysis",
    theme = theme (plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
    plot.margin = margin(t =10, r = 10, b = 10, l = 10)))


ggsave("Desktop/full.pdf", combined_plot, width = 28, height = 24)


