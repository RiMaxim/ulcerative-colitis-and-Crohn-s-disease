# ==============================================================================
# УПРОЩЕННЫЙ СКРИПТ: КОРРЕЛЯЦИОННЫЙ АНАЛИЗ СО ВСЕМИ МЕТАДАННЫМИ
# ==============================================================================

library(dplyr)
library(reshape2)
library(ggplot2)
library(tidyr)

# ==============================================================================
# 1. ЗАГРУЗКА ДАННЫХ
# ==============================================================================

data_rel <- read.table("Desktop/WORK/gut/1_stage/R/input.tsv", 
                       header = TRUE, sep = "\t", 
                       check.names = FALSE, fill = TRUE)

metadata <- read.table("Desktop/WORK/gut/1_stage/R/metadata.tsv", 
                       header = TRUE, sep = "\t", 
                       check.names = FALSE, fill = TRUE)

cat("=== ЗАГРУЖЕННЫЕ ДАННЫЕ ===\n")
cat("Размер data_rel:", dim(data_rel), "\n")
cat("Размер metadata:", dim(metadata), "\n")
str(metadata)
# ==============================================================================
# 2. ВЫРАВНИВАНИЕ ОБРАЗЦОВ И ФИЛЬТРАЦИЯ IBD
# ==============================================================================

# Выравнивание по Sample ID
ids_rel  <- as.character(data_rel[, 1])
ids_meta <- as.character(metadata[, 1])
common_samples <- intersect(ids_meta, ids_rel)

if (length(common_samples) == 0) {
  stop("Ошибка: совпадений по Sample ID не найдено!")
}

metadata_aligned <- metadata[match(common_samples, ids_meta), , drop = FALSE]
data_rel_aligned <- data_rel[match(common_samples, ids_rel), , drop = FALSE]

# Матрица относительных обилий
mat_rel <- as.matrix(sapply(data_rel_aligned[, -1], as.numeric))
rownames(mat_rel) <- common_samples

# Исправляем безымянную первую колонку
colnames(metadata_aligned)[1] <- "sample_id"

# Фильтруем только IBD пациентов (CD и UC)
meta_ibd <- metadata_aligned %>%
  filter(diagnosis %in% c("CD", "UC"))

# Фильтруем матрицу бактерий
data_ibd <- mat_rel[rownames(mat_rel) %in% meta_ibd$sample_id, ]

cat("\n=== ФИЛЬТРАЦИЯ ===\n")
cat(sprintf("IBD пациентов: %d\n", nrow(meta_ibd)))
cat(sprintf("  CD: %d\n", sum(meta_ibd$diagnosis == "CD")))
cat(sprintf("  UC: %d\n", sum(meta_ibd$diagnosis == "UC")))

# ==============================================================================
# 3. ОТБОР ТОП-20 ВИДОВ
# ==============================================================================

n_species <- ncol(data_ibd)
k_top <- min(1000, n_species)
top_species <- names(sort(colMeans(data_ibd), decreasing = TRUE)[1:k_top])
data_ibd_top <- data_ibd[, top_species, drop = FALSE]

cat("\n=== ТОП-20 ВИДОВ ===\n")
print(top_species)

# ==============================================================================
# 4. ОБРАБОТКА ТЕРАПИИ (комбинированные показатели)
# ==============================================================================

drug_cols <- grep("^drug_", colnames(meta_ibd), value = TRUE)

if (length(drug_cols) > 0) {
  meta_ibd <- meta_ibd %>%
    mutate(across(all_of(drug_cols), ~ as.numeric(. == "1"))) %>%
    mutate(
      therapy_count = rowSums(select(., all_of(drug_cols)), na.rm = TRUE),
      any_therapy = ifelse(therapy_count > 0, 1, 0)
    )
} else {
  # Если колонок с лекарствами нет, создаем пустые
  meta_ibd$therapy_count <- NA
  meta_ibd$any_therapy <- NA
}

# ==============================================================================
# 5. ПОЛНЫЙ СПИСОК КЛИНИЧЕСКИХ ПАРАМЕТРОВ
# ==============================================================================

clinical_params <- c(
  # Демография
  "age", "sex",
  
  # Воспалительные маркеры
  "fc", "crp", "lbp", "wbc", "hb", "rbc",
  
  # Аутоантитела
  "asca_igg", "asca_iga", "anca_igg", "anca_iga",
  
  # SCFA абсолютные
  "c2_abs", "c3_abs", "c4_abs", "sum_isocn_abs",
  
  # SCFA относительные
  "c2", "c3", "c4", "sum_isocn",
  
  # Метаболические индексы
  "anaer_index", "bile_acid", "isocn_cn", "isoc5_c5", "sum",
  
  # Эндоскопия - терминальный подвздошная кишка
  "terminal_ileum_1", "terminal_ileum_2", "terminal_ileum_3",
  "terminal_ileum_4", "terminal_ileum_5", "terminal_ileum_6",
  "terminal_ileum_7", "terminal_ileum_8", "terminal_ileum_9",
  "terminal_ileum_10", "terminal_ileum_11",
  
  # Эндоскопия - слепая кишка
  "cecal_pole_1", "cecal_pole_2", "cecal_pole_3", "cecal_pole_4",
  "cecal_pole_5", "cecal_pole_6", "cecal_pole_7", "cecal_pole_8",
  "cecal_pole_9", "cecal_pole_10", "cecal_pole_11",
  
  # Терапия
  "therapy_count", "any_therapy"
)

# Проверяем, какие параметры реально существуют в данных
existing_params <- clinical_params[clinical_params %in% colnames(meta_ibd)]

cat("\n=== ДОСТУПНЫЕ КЛИНИЧЕСКИЕ ПАРАМЕТРЫ ===\n")
cat(sprintf("Всего параметров: %d\n", length(clinical_params)))
cat(sprintf("Доступно в данных: %d\n", length(existing_params)))

# Показываем отсутствующие параметры
missing_params <- setdiff(clinical_params, existing_params)
if (length(missing_params) > 0) {
  cat("\nОтсутствующие параметры:\n")
  print(missing_params)
}

# ==============================================================================
# 6. ФУНКЦИЯ ПРЕОБРАЗОВАНИЯ ТИТРОВ
# ==============================================================================

convert_titer <- function(x) {
  if (is.na(x) || x %in% c("ND", "", "ND ", "NA", "none", "N/A")) return(NA)
  if (grepl("^\\d+:\\d+$", x)) {
    return(as.numeric(strsplit(x, ":")[[1]][2]))
  }
  if (grepl("^\\d+/\\d+$", x)) {
    return(as.numeric(strsplit(x, "/")[[1]][2]))
  }
  return(suppressWarnings(as.numeric(x)))
}

# ==============================================================================
# 7. ПОДГОТОВКА КЛИНИЧЕСКИХ ДАННЫХ
# ==============================================================================

cat("\n=== ПОДГОТОВКА КЛИНИЧЕСКИХ ДАННЫХ ===\n")

# Определяем категории колонок для специальной обработки
autoantibody_cols <- intersect(c("asca_igg", "asca_iga", "anca_igg", "anca_iga"), existing_params)
endoscopy_cols <- intersect(grep("terminal_ileum_|cecal_pole_", existing_params, value = TRUE), existing_params)

# Преобразуем данные
meta_clean <- meta_ibd %>%
  select(all_of(existing_params)) %>%
  mutate(
    # Аутоантитела - преобразование титров
    across(all_of(autoantibody_cols), ~ sapply(as.character(.), convert_titer)),
    
    # Пол: M=1, F=0
    sex = ifelse(sex == "M", 1, ifelse(sex == "F", 0, NA)),
    
    # Эндоскопия - оставляем как есть (уже числа 0,1,2,3)
    across(all_of(endoscopy_cols), ~ suppressWarnings(as.numeric(as.character(.)))),
    
    # Все остальные параметры - просто числа
    across(-all_of(c(autoantibody_cols, "sex", endoscopy_cols)), 
           ~ suppressWarnings(as.numeric(as.character(.))))
  )

cat(sprintf("Размер очищенных данных: %d x %d\n", nrow(meta_clean), ncol(meta_clean)))

# Статистика по пропускам
na_counts <- sapply(meta_clean, function(x) sum(is.na(x)))
cat("\nПараметры с пропусками (>10%):\n")
print(names(na_counts[na_counts > nrow(meta_clean) * 0.1]))

# ==============================================================================
# 8. РАСЧЕТ КОРРЕЛЯЦИЙ (С ЗАЩИТОЙ ОТ ЛОЖНЫХ КОРРЕЛЯЦИЙ) 
# ==============================================================================

cat("\n=== РАСЧЕТ КОРРЕЛЯЦИЙ ===\n")
cat(sprintf("Видов бактерий: %d\n", length(top_species)))
cat(sprintf("Клинических параметров: %d\n", ncol(meta_clean)))
cat(sprintf("Всего пар: %d\n", length(top_species) * ncol(meta_clean)))

# Матрицы для результатов
cor_matrix <- matrix(
  NA_real_,
  nrow = length(top_species),
  ncol = ncol(meta_clean),
  dimnames = list(top_species, colnames(meta_clean))
)

p_matrix <- matrix(
  NA_real_,
  nrow = length(top_species),
  ncol = ncol(meta_clean),
  dimnames = list(top_species, colnames(meta_clean))
)

method_matrix <- matrix(
  NA_character_,
  nrow = length(top_species),
  ncol = ncol(meta_clean),
  dimnames = list(top_species, colnames(meta_clean))
)


# Счетчики
filtered_pairs <- 0
total_pairs <- length(top_species) * ncol(meta_clean)
pair_count <- 0

for (sp in top_species) {
  for (cl in colnames(meta_clean)) {
    pair_count <- pair_count + 1
    
    vec_sp <- data_ibd_top[, sp]
    vec_cl <- meta_clean[, cl]
    
    # 1. Удаление NA пар
    complete <- complete.cases(vec_sp, vec_cl)
    n_complete <- sum(complete)
    
    # Минимальный порог наблюдений для запуска теста
    if (n_complete < 10) {
      filtered_pairs <- filtered_pairs + 1
      next
    }
    
    vec_sp_clean <- vec_sp[complete]
    vec_cl_clean <- vec_cl[complete]
    
    # 2. Фильтр редкости бактерий (обилие > 0)
    # Порог изменен до 6-8 присутствий, чтобы не отсечь важные редкие виды
    n_present_sp <- sum(vec_sp_clean > 0)
    if (n_present_sp < 7) { 
      filtered_pairs <- filtered_pairs + 1
      next
    }
    
    # 3. Проверка на вариабельность признаков (константность)
    unique_cl <- unique(vec_cl_clean)
    unique_sp <- unique(vec_sp_clean)
    
    if (length(unique_sp) > 1 && length(unique_cl) > 1) {
      
      # ОПРЕДЕЛЕНИЕ МЕТОДА: бинарный или многоуровневый
      # Если у параметра ровно 2 исхода (например, 0 и 1), используем Пирсона (точечно-бисериальную)
      if (length(unique_cl) == 2) {
        current_method <- "pearson"
      } else {
        current_method <- "spearman"
      }
      
      tryCatch({
        if (current_method == "pearson") {
          # Точечно-бисериальная корреляция
          test <- cor.test(vec_sp_clean, as.numeric(vec_cl_clean), method = "pearson")
        } else {
          # Ранговая корреляция Спирмена с exact=FALSE для обработки связанных рангов (ties)
          test <- cor.test(vec_sp_clean, vec_cl_clean, method = "spearman", exact = FALSE)
        }
        
        cor_matrix[sp, cl] <- test$estimate
        p_matrix[sp, cl] <- test$p.value
        method_matrix[sp, cl] <- current_method
        
      }, error = function(e) {
        cor_matrix[sp, cl] <- NA
        p_matrix[sp, cl] <- NA
      })
      
    } else {
      filtered_pairs <- filtered_pairs + 1
    }
    
    if (pair_count %% 100 == 0) {
      cat(sprintf("Прогресс: %.1f%% (%d/%d), отфильтровано: %d\n", 
                  pair_count/total_pairs*100, pair_count, total_pairs, filtered_pairs))
    }
  }
}

cat("Расчет матриц завершен!\n")
cat(sprintf("Отфильтровано пар: %d\n", filtered_pairs))

# 4. Поправка на множественное тестирование (FDR)
cat("Применение поправки Бенджамини-Хохберга (FDR)...\n")
p_adj_vector <- p.adjust(as.vector(p_matrix), method = "BH")
p_adj_matrix <- matrix(p_adj_vector, 
                       nrow = nrow(p_matrix), 
                       ncol = ncol(p_matrix), 
                       dimnames = dimnames(p_matrix))
df_padj <- melt(
  p_adj_matrix,
  varnames = c("Species", "Parameter"),
  value.name = "Padj"
)

# ==============================================================================
# 9. ПОДГОТОВКА ДЛЯ ВИЗУАЛИЗАЦИИ (С FDR КОРРЕКЦИЕЙ)
# ==============================================================================

# Преобразуем в длинный формат
df_cor <- melt(cor_matrix, varnames = c("Species", "Parameter"), value.name = "Correlation")
df_p <- melt(p_matrix, varnames = c("Species", "Parameter"), value.name = "Pvalue")
df_plot <- merge(
  df_cor,
  df_p,
  by = c("Species", "Parameter")
) %>%
  left_join(
    df_padj,
    by = c("Species", "Parameter")
  ) %>%
  mutate(
    
    # Звездочки для p-value и p.adj раздельно
    Stars_pvalue = case_when(
      Pvalue < 0.001 ~ "***",
      Pvalue < 0.01 ~ "**",
      Pvalue < 0.05 ~ "*",
      TRUE ~ ""
    ),
    
    Stars_padj = case_when(
      Padj < 0.001 ~ "***",
      Padj < 0.01 ~ "**",
      Padj < 0.05 ~ "*",
      TRUE ~ ""
    ),
    
    # Комбинируем звездочки: p-value / p.adj
    Significance = ifelse(
      Stars_pvalue != "" | Stars_padj != "",
      paste(Stars_pvalue, Stars_padj, sep = "/"),
      ""
    ),
    
    # Остальные категории как были
    Category = case_when(
      Parameter %in% c("age", "sex") ~ "1_Demographics",
      Parameter %in% c("fc", "crp", "lbp", "wbc", "hb", "rbc") ~ "2_Inflammation",
      Parameter %in% c("asca_igg", "asca_iga", "anca_igg", "anca_iga") ~ "3_Serology",
      Parameter %in% c("c2_abs", "c3_abs", "c4_abs", "sum_isocn_abs") ~ "4_SCFA absolute values",
      Parameter %in% c("c2", "c3", "c4", "sum_isocn") ~ "4_SCFA relative values",
      Parameter %in% c("isocn_cn", "isoc5_c5", "sum", "anaer_index") ~ "5_Ratios & Indices",
      Parameter %in% c("bile_acid") ~ "6_FBA",
      grepl("terminal_ileum_", Parameter) ~ "7_Terminal Ileum",
      grepl("cecal_pole_", Parameter) ~ "8_Cecal Pole",
      Parameter %in% c("therapy_count", "any_therapy") ~ "Therapy",
      TRUE ~ "10_Other"
    )
  )


# Удаляем NA
df_plot <- df_plot %>% filter(!is.na(Correlation))

cat(sprintf("\nПодготовлено %d пар для визуализации\n", nrow(df_plot)))
cat(sprintf("Значимых корреляций (FDR < 0.05): %d\n", sum(df_plot$Padj < 0.05, na.rm = TRUE)))

# ==============================================================================
# 9.4. ФИЛЬТРАЦИЯ: ТОЛЬКО ВИДЫ С ХОТЯ БЫ ОДНИМ * 
# ==============================================================================

cat("\n=== ФИЛЬТРАЦИЯ ВИДОВ СО ЗНАЧИМЫМИ КОРРЕЛЯЦИЯМИ ===\n")

# Находим виды, у которых есть хотя бы одна звездочка
species_with_stars <- df_plot %>%
  filter(Pvalue < 0.05) %>% #Padj
  pull(Species) %>%
  unique()

cat(sprintf("Видов до фильтрации: %d\n", length(unique(df_plot$Species))))
cat(sprintf("Видов с Pvalue < 0.05: %d\n", #FDR
            length(species_with_stars)))

# Фильтруем df_plot
df_plot_filtered <- df_plot %>%
  filter(Species %in% species_with_stars)

cat(sprintf("Строк после фильтрации: %d\n", nrow(df_plot_filtered)))

# Обновляем df_plot для дальнейшего использования
df_plot <- df_plot_filtered

# Если нужно также обновить cor_matrix для кластеризации
if (length(species_with_stars) > 1) {
  cor_matrix <- cor_matrix[species_with_stars, , drop = FALSE]
  cat("Матрица корреляций обновлена для кластеризации\n")
} else if (length(species_with_stars) == 1) {
  cat("ВНИМАНИЕ: Только 1 вид со звездочками, кластеризация не возможна\n")
} else {
  cat("ВНИМАНИЕ: Нет видов со звездочками, отображение невозможно\n")
}

# ==============================================================================
# 9.5. КЛАСТЕРИЗАЦИЯ ПО ОСИ Y (ВИДЫ БАКТЕРИЙ)
# ==============================================================================

cat("\n=== КЛАСТЕРИЗАЦИЯ ВИДОВ ===\n")

# Создаем матрицу корреляций для кластеризации видов
# Используем средние корреляции для каждого вида
species_cor_matrix <- cor_matrix

# Расчет расстояния между видами на основе их корреляций со всеми параметрами
# Заменяем NA на 0 для кластеризации
species_dist_matrix <- species_cor_matrix
species_dist_matrix[is.na(species_dist_matrix)] <- 0

# Вычисляем евклидово расстояние между видами
if (nrow(species_dist_matrix) < 2) {
  
  stop(
    "После фильтрации осталось менее двух видов. Кластеризация невозможна."
  )
  
}
species_dist <- dist(species_dist_matrix, method = "euclidean")
species_clust <- hclust(species_dist, method = "ward.D2")

# Получаем порядок видов для кластеризации
species_order <- species_clust$labels[species_clust$order]

cat("Порядок видов после кластеризации:\n")
print(species_order)

# Преобразуем Species в фактор с порядком кластеризации
df_plot$Species <- factor(df_plot$Species, levels = species_order)

# ==============================================================================
# 10. ХИТМЭП (С КЛАСТЕРИЗАЦИЕЙ ПО ОСИ Y)
# ==============================================================================

# Упорядочиваем категории
df_plot$Category <- factor(df_plot$Category, 
                           levels = sort(unique(df_plot$Category)))

# Ширина хитмэпа зависит от количества параметров
plot_width <- max(12, ncol(meta_clean) * 0.3)

# В секции 10, измените subtitle в p_heatmap:

p_heatmap <- ggplot(df_plot, aes(x = Parameter, y = Species, fill = Correlation)) +
  geom_tile(color = "white", linewidth = 0.3) +
  geom_text(aes(label = Significance), color = "black", size = 2.5, vjust = 0.7) +
  scale_fill_gradient2(low = "#4575b4", mid = "#ffffbf", high = "#d73027", 
                       midpoint = 0, limits = c(-1, 1),
                       name = expression(rho)) +
  facet_grid(. ~ Category, scales = "free_x", space = "free") +
  labs(
    title = "Gut Microbiota Correlations with All Clinical Parameters in IBD",
    subtitle = sprintf(
      "Pearson correlations for binary variables (0/1) and Spearman correlations for continuous or ordinal variables. IBD patients (n=%d). Significance: * < 0.05, ** < 0.01, *** < 0.001. Format: p-value / FDR-adjusted p-value",
      nrow(data_ibd_top)
    ),
    x = "Clinical Parameters",
    y = "Bacterial Species"
  ) +
  theme_minimal(base_size = 9) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    axis.text.y = element_text(face = "italic", size = 12),
    axis.title.x   = element_text(size = 14),
    axis.title.y   = element_text(size = 14),
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    plot.subtitle = element_text(size = 9, hjust = 0.5),
    strip.text = element_text(face = "bold", size = 8),
    strip.background = element_rect(fill = "lightgray", color = NA),
    panel.grid = element_blank()
  )

# Сохраняем
pdf("Desktop/WORK/gut/1_stage/R/correlations_all_params.pdf", 
    width = 20, height = 8)
print(p_heatmap)
dev.off()




###################################Проверка

sp <- "Phocaeicola_coprocola"
cl <- "cecal_pole_9"

# берём ровно те же данные, что в цикле
vec_sp <- data_ibd_top[, sp]
vec_cl <- meta_clean[, cl]

# применяем ТО ЖЕ условие фильтрации
complete <- complete.cases(vec_sp, vec_cl)

# итоговые данные, которые реально пошли в cor.test
cor_data <- data.frame(
  sample_id = rownames(data_ibd_top)[complete],
  Phocaeicola_coprocola = vec_sp[complete],
  cecal_pole_9 = vec_cl[complete]
)

cor_data

cor.test(
  cor_data$Phocaeicola_coprocola,
  cor_data$cecal_pole_9,
  method = "pearson",
  exact = FALSE
)


################################################################################Таблица для графика



df_plot[abs(df_plot$Padj) < 0.05, ]


