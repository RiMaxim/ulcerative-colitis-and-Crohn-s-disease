library(tidyverse)
library(corrplot)
library(pheatmap)
library(ggplot2)
library(ggpubr)
library(Hmisc)
library(psych)
library(vegan)
library(mixOmics)
library(ComplexHeatmap)
library(compositions)
library(caret) 
library(randomForest)
library(ROCR)
library(reshape2)
library(reshape2)

input_data <- read.table("Desktop/WORK/gut/1_stage/R/input.tsv", header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)
input_metadata <- read.table("Desktop/WORK/gut/1_stage/R/metadata.tsv", header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)
input_fecal <- read.table("Desktop/WORK/gut/1_stage/R/fecal.txt", header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)
input_lipidome <- read.table("Desktop/WORK/gut/1_stage/R/lipidome.txt", header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)
input_serum <- read.table("Desktop/WORK/gut/1_stage/R/serum.txt", header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)


diagnosis_col  <- "diagnosis"



# --- Шаг 3: Функции для работы с данными ---

get_measured_samples <- function(df) {
  if(ncol(df) == 0) return(character(0))
  first_col <- trimws(as.character(df[, 1]))
  rownames(df)[first_col != "ND" & 
                 !is.na(first_col) & 
                 first_col != "" & 
                 first_col != "NA" & 
                 first_col != "NULL"]
}

to_numeric_matrix <- function(df) {
  df[df == "ND" | df == "" | df == "NA" | df == "NULL"] <- NA
  mat <- matrix(as.numeric(as.matrix(df)), nrow = nrow(df), 
                dimnames = list(rownames(df), colnames(df)))
  for(col in 1:ncol(mat)) {
    if(any(is.na(mat[, col]))) {
      med_val <- median(mat[, col], na.rm = TRUE)
      if(is.na(med_val)) med_val <- 0
      mat[is.na(mat[, col]), col] <- med_val
    }
  }
  return(mat)
}

preprocess_scaled <- function(mat) {
  mat_log <- log2(mat + 1)
  mat_log[is.infinite(mat_log) | is.nan(mat_log) | is.na(mat_log)] <- 0
  return(scale(mat_log))
}

remove_zero_var <- function(mat, block_name) {
  if(is.null(mat) || ncol(mat) == 0) return(mat)
  if(any(is.na(mat))) mat[is.na(mat)] <- 0
  if(any(is.infinite(mat))) mat[is.infinite(mat)] <- 0
  
  variances <- apply(mat, 2, function(x) var(x, na.rm = TRUE))
  if(any(is.na(variances))) variances[is.na(variances)] <- 0
  
  zero_var <- variances == 0
  if(sum(zero_var, na.rm = TRUE) > 0) {
    cat("Блок", block_name, "- удалено признаков:", sum(zero_var, na.rm = TRUE), "\n")
    return(mat[, !zero_var, drop = FALSE])
  } else {
    cat("Блок", block_name, "- признаков с нулевой дисперсией не найдено\n")
    return(mat)
  }
}

# --- Шаг 4: Подготовка данных ---

cat("\n--- Подготовка микробиома ---\n")

# Микробиом
samples_micro <- rownames(input_metadata)[!is.na(input_metadata[[diagnosis_col]]) & 
                                            input_metadata[[diagnosis_col]] != "ND" &
                                            input_metadata[[diagnosis_col]] != ""]
samples_micro <- intersect(samples_micro, rownames(input_data))

meta_micro <- input_metadata[samples_micro, ]
micro_samples <- rownames(meta_micro)[meta_micro[[diagnosis_col]] %in% c("CD", "UC")]

cat("Образцов CD/UC для микробиома:", length(micro_samples), "\n")

data_micro <- input_data[micro_samples, ]
meta_micro_clean <- input_metadata[micro_samples, ]

# CLR трансформация микробиома
keep_taxa <- colSums(data_micro > 0) > (0.10 * nrow(data_micro))
micro_clr <- as.matrix(clr(data_micro[, keep_taxa] + 1))
micro_clr[is.infinite(micro_clr) | is.nan(micro_clr)] <- 0
micro_clr <- remove_zero_var(micro_clr, "microbiome")

cat("Размер микробиома:", dim(micro_clr), "\n")

cat("\n--- Подготовка фекальных метаболитов ---\n")




# Фекальные метаболиты
samples_fecal <- get_measured_samples(input_fecal)
samples_fecal <- intersect(samples_fecal, rownames(input_metadata))
meta_fecal <- input_metadata[samples_fecal, ]
fecal_samples <- rownames(meta_fecal)[meta_fecal[[diagnosis_col]] %in% c("CD", "UC")]

cat("Образцов CD/UC для фекальных метаболитов:", length(fecal_samples), "\n")

fecal_num <- to_numeric_matrix(input_fecal[fecal_samples, ])
fecal_ready <- preprocess_scaled(fecal_num)
fecal_ready <- remove_zero_var(fecal_ready, "fecal")

cat("Размер фекальных:", dim(fecal_ready), "\n")

# --- Шаг 5: Находим общие образцы ---

common_samples <- intersect(micro_samples, fecal_samples)
cat("\nОбщих образцов:", length(common_samples), "\n")

# Фильтруем данные по общим образцам
micro_common <- micro_clr[common_samples, ]
fecal_common <- fecal_ready[common_samples, ]
meta_common <- input_metadata[common_samples, ]

# Создаем фактор диагноза
y_common <- factor(meta_common[[diagnosis_col]])
y_common <- droplevels(y_common)

cat("\nРаспределение классов:\n")
print(table(y_common))

# --- Шаг 6: Разделение на CD и UC ---

cd_idx <- which(y_common == "CD")
uc_idx <- which(y_common == "UC")

cat("\nCD образцов:", length(cd_idx), "\n")
cat("UC образцов:", length(uc_idx), "\n")

# Данные для CD
micro_cd <- micro_common[cd_idx, ]
fecal_cd <- fecal_common[cd_idx, ]

# Данные для UC
micro_uc <- micro_common[uc_idx, ]
fecal_uc <- fecal_common[uc_idx, ]

# --- Шаг 7: ИСПРАВЛЕННАЯ ФУНКЦИЯ для расчета корреляций ---


# ============================================================================
# ПРОСТАЯ ФУНКЦИЯ КОРРЕЛЯЦИИ СПИРМЕНА
# ============================================================================

calculate_spearman <- function(micro_data, fecal_data, group_name) {
  
  cat("\n=========================================")
  cat("\nКорреляции Спирмена для:", group_name)
  cat("\n=========================================\n")
  
  n_bacteria <- ncol(micro_data)
  n_metabolites <- ncol(fecal_data)
  
  cat("  Бактерий:", n_bacteria, "\n")
  cat("  Метаболитов:", n_metabolites, "\n")
  cat("  Всего пар:", n_bacteria * n_metabolites, "\n")
  
  # Матрицы для результатов
  cor_matrix <- matrix(NA, nrow = n_bacteria, ncol = n_metabolites)
  p_matrix <- matrix(NA, nrow = n_bacteria, ncol = n_metabolites)
  
  rownames(cor_matrix) <- colnames(micro_data)
  colnames(cor_matrix) <- colnames(fecal_data)
  rownames(p_matrix) <- colnames(micro_data)
  colnames(p_matrix) <- colnames(fecal_data)
  
  # Считаем корреляции
  cat("\n  Расчет...\n")
  
  for(i in 1:n_bacteria) {
    for(j in 1:n_metabolites) {
      # Просто считаем корреляцию
      test <- cor.test(micro_data[, i], fecal_data[, j], method = "spearman")
      
      cor_matrix[i, j] <- test$estimate
      p_matrix[i, j] <- test$p.value
    }
    
    # Прогресс
    if(i %% 10 == 0) cat("    Обработано", i, "из", n_bacteria, "бактерий\n")
  }
  
  # Создаем таблицу результатов
  results <- data.frame(
    Bacteria = rep(rownames(cor_matrix), each = n_metabolites),
    Metabolite = rep(colnames(cor_matrix), times = n_bacteria),
    Correlation = as.vector(cor_matrix),
    P_value = as.vector(p_matrix),
    Group = group_name
  )
  
  # Убираем NA
  results <- results[!is.na(results$Correlation), ]
  
  # Флаг значимости
  results$Significant <- results$P_value < 0.05
  
  # Сортируем
  results <- results[order(-abs(results$Correlation)), ]
  
  # Топ-10
  cat("\n  ТОП-10 КОРРЕЛЯЦИЙ:\n")
  print(head(results[, c("Bacteria", "Metabolite", "Correlation", "P_value")], 10))
  
  cat("\n  Значимых (p < 0.05):", sum(results$Significant), "из", nrow(results), "\n")
  
  return(list(
    cor_matrix = cor_matrix,
    p_matrix = p_matrix,
    results = results
  ))
}

# ============================================================================
# ЗАПУСК
# ============================================================================

cat("\n=========================================")
cat("\nЗАПУСК АНАЛИЗА")
cat("\n=========================================\n")

# Для CD
res_cd <- calculate_spearman(micro_cd, fecal_cd, "CD")

# Для UC
res_uc <- calculate_spearman(micro_uc, fecal_uc, "UC")

# Для всех
res_all <- calculate_spearman(micro_common, fecal_common, "ALL")

# ============================================================================
# СОХРАНЕНИЕ
# ============================================================================

cat("\n=========================================")
cat("\nСОХРАНЕНИЕ РЕЗУЛЬТАТОВ")
cat("\n=========================================\n")

# Сохраняем результаты
write.table(res_cd$results, "correlations_CD.tsv", sep = "\t", row.names = FALSE)
write.table(res_uc$results, "correlations_UC.tsv", sep = "\t", row.names = FALSE)
write.table(res_all$results, "correlations_ALL.tsv", sep = "\t", row.names = FALSE)

cat("  correlations_CD.tsv\n")
cat("  correlations_UC.tsv\n")
cat("  correlations_ALL.tsv\n")

# Только значимые
write.table(res_cd$results[res_cd$results$Significant, ], 
            "correlations_significant_CD.tsv", sep = "\t", row.names = FALSE)
write.table(res_uc$results[res_uc$results$Significant, ], 
            "correlations_significant_UC.tsv", sep = "\t", row.names = FALSE)

cat("  correlations_significant_CD.tsv\n")
cat("  correlations_significant_UC.tsv\n")




# ============================================================================
# ТЕПЛОВЫЕ КАРТЫ
# ============================================================================

# --- ФУНКЦИЯ ДЛЯ СОЗДАНИЯ ТЕПЛОВОЙ КАРТЫ ---

create_heatmap_ggplot <- function(cor_matrix, title, filename) {
  
  # Преобразуем в long format
  cor_melt <- melt(cor_matrix, na.rm = FALSE)  # na.rm = FALSE чтобы сохранить NA
  colnames(cor_melt) <- c("Bacteria", "Metabolite", "Correlation")
  
  # Заменяем _ на пробел
  cor_melt$Bacteria <- gsub("_", " ", cor_melt$Bacteria)
  cor_melt$Metabolite <- gsub("_", " ", cor_melt$Metabolite)
  
  # Создаем тепловую карту
  p <- ggplot(cor_melt, aes(x = Metabolite, y = Bacteria, fill = Correlation)) +
    # ★★★ БЕЛЫЕ ПОЛОСКИ МЕЖДУ ЯЧЕЙКАМИ ★★★
    geom_tile(color = "white", linewidth = 0.5) +  
    # ★★★ СЕРЫЙ ЦВЕТ ДЛЯ NA ★★★
    scale_fill_gradient2(
      low = "darkblue", 
      mid = "white", 
      high = "darkred", 
      midpoint = 0, 
      limits = c(-1, 1),
      na.value = "gray"  # ← СЕРЫЙ ДЛЯ NA
    ) +
    geom_text(aes(label = ifelse(is.na(Correlation), "", sprintf("%.2f", Correlation))), 
              size = 2, na.rm = TRUE) +
    theme_minimal() +
    theme(
      axis.text.y = element_text(face = "italic", size = 10),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
      axis.title = element_blank(),
      panel.grid = element_blank(),
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      # ★★★ УБИРАЕМ ВСЕ ЛИШНИЕ ЛИНИИ ★★★
      panel.border = element_blank(),
      axis.ticks = element_blank()
    ) +
    labs(title = title)
  
  # Сохраняем
  ggsave(filename, p, width = 8, height = 5.5)
  cat("  ✅", filename, "\n")
}

# ============================================================================
# РИСУНОК 1: CD
# ============================================================================

cat("\n--- Рисунок 1: heatmap_CD_strong.pdf ---\n")

cor_cd <- res_cd$cor_matrix
p_cd <- res_cd$p_matrix

# Фильтрация
cor_cd[p_cd >= 0.05] <- NA
cor_cd[abs(cor_cd) <= 0.5] <- NA

# Удаляем пустые строки и столбцы
cor_cd <- cor_cd[rowSums(!is.na(cor_cd)) > 0, , drop = FALSE]
cor_cd <- cor_cd[, colSums(!is.na(cor_cd)) > 0, drop = FALSE]

cat("  Осталось бактерий:", nrow(cor_cd), "\n")
cat("  Осталось метаболитов:", ncol(cor_cd), "\n")
cat("  Всего корреляций:", sum(!is.na(cor_cd)), "\n")

if(nrow(cor_cd) > 0 && ncol(cor_cd) > 0) {
  create_heatmap_ggplot(cor_cd, 
                        "Fecal metabolome VS. Microbiome\nCD, |Rho| > 0.5, P_value < 0.05", 
                        "Desktop/WORK/gut/1_stage/R/heatmap_CD_strong.pdf")
} else {
  cat("  ⚠️ Нет данных для CD\n")
}


# ============================================================================
# РИСУНОК 2: UC
# ============================================================================

cat("\n--- Рисунок 2: heatmap_UC_strong.pdf ---\n")

cor_uc <- res_uc$cor_matrix
p_uc <- res_uc$p_matrix

cor_uc[p_uc >= 0.05] <- NA
cor_uc[abs(cor_uc) <= 0.6] <- NA

cor_uc <- cor_uc[rowSums(!is.na(cor_uc)) > 0, , drop = FALSE]
cor_uc <- cor_uc[, colSums(!is.na(cor_uc)) > 0, drop = FALSE]

cat("  Осталось бактерий:", nrow(cor_uc), "\n")
cat("  Осталось метаболитов:", ncol(cor_uc), "\n")
cat("  Всего корреляций:", sum(!is.na(cor_uc)), "\n")

if(nrow(cor_uc) > 0 && ncol(cor_uc) > 0) {
  create_heatmap_ggplot(cor_uc, 
                        "Fecal metabolome VS. Microbiome\nUC, |Rho| > 0.6, P_value < 0.05", 
                        "Desktop/WORK/gut/1_stage/R/heatmap_UC_strong.pdf")
} else {
  cat("  ⚠️ Нет данных для UC\n")
}
















































##################beta





# --- Шаг 1: Загрузка библиотек и данных ---
library(vegan)         # Для adonis2, ordiellipse
library(compositions)  # Для функции clr()

input_data <- read.table("Desktop/WORK/gut/1_stage/R/input.tsv", header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)
input_metadata <- read.table("Desktop/WORK/gut/1_stage/R/metadata.tsv", header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)
input_fecal <- read.table("Desktop/WORK/gut/1_stage/R/fecal.txt", header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)
input_lipidome <- read.table("Desktop/WORK/gut/1_stage/R/lipidome.txt", header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)
input_serum <- read.table("Desktop/WORK/gut/1_stage/R/serum.txt", header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)

diagnosis_col  <- "diagnosis"

# --- Шаг 2: Синхронизация и фильтрация образцов ---
cat("\n--- 2. Фильтрация структурных пропусков ---\n")

get_measured_samples <- function(df) {
  rownames(df)[df[, 1] != "ND" & !is.na(df[, 1]) & df[, 1] != ""]
}

samples_fecal <- get_measured_samples(input_fecal)
samples_lipid <- get_measured_samples(input_lipidome)
samples_serum <- get_measured_samples(input_serum)
samples_diag  <- rownames(input_metadata)[!is.na(input_metadata[[diagnosis_col]]) & input_metadata[[diagnosis_col]] != "ND"]


complete_samples <- intersect(samples_diag, samples_fecal)
complete_samples <- intersect(samples_diag, samples_lipid)
complete_samples <- intersect(samples_diag, samples_serum)

metadata_sub <- input_metadata[complete_samples, ]
ibd_samples  <- rownames(metadata_sub)[metadata_sub[[diagnosis_col]] %in% c("CD", "UC")]

meta_clean  <- input_metadata[ibd_samples, ]
data_clean  <- input_data[ibd_samples, ]
fecal_clean <- input_fecal[ibd_samples, ]
lipid_clean <- input_lipidome[ibd_samples, ]
serum_clean <- input_lipidome[ibd_samples, ]

# Переименование для caret
y_raw <- factor(meta_clean[[diagnosis_col]])
y_raw <- droplevels(y_raw)
levels(y_raw) <- c("Class_CD", "Class_UC")
y_model <- factor(y_raw, levels = c("Class_CD", "Class_UC"))

cat("Итоговая выборка (CD и UC):", length(y_model), "\n")
print(table(y_model))

# --- Шаг 3: Конвертация типов и нормализация данных ---
cat("\n--- 3. Биологическая нормализация и предобработка ---\n")

to_numeric_matrix <- function(df) {
  df[df == "ND"] <- NA
  mat <- matrix(as.numeric(as.matrix(df)), nrow = nrow(df), dimnames = dimnames(df))
  for(col in 1:ncol(mat)) {
    if(any(is.na(mat[, col]))) mat[is.na(mat[, col]), col] <- median(mat[, col], na.rm = TRUE)
  }
  return(mat)
}

fecal_num <- to_numeric_matrix(fecal_clean)
lipid_num <- to_numeric_matrix(lipid_clean)
serum_num <- to_numeric_matrix(serum_clean)

# Микробиом: Очистка шума (10%) + CLR-трансформация
keep_taxa <- colSums(data_clean > 0) > (0.10 * nrow(data_clean))
# Преобразование в матрицу для корректной работы clr
micro_clr <- as.matrix(clr(as.matrix(data_clean[, keep_taxa]) + 1))

# Вещества: Log2 + Z-масштабирование
preprocess_scaled <- function(mat) {
  mat_log <- log2(mat + 1)
  mat_log[is.infinite(mat_log) | is.nan(mat_log) | is.na(mat_log)] <- 0
  return(scale(mat_log))
}

fecal_ready <- preprocess_scaled(fecal_num)
lipid_ready <- preprocess_scaled(lipid_num)
serum_ready <- preprocess_scaled(serum_num)

# --- Шаг 4: Анализ бета-разнообразия фекального метаболома ---
cat("\n--- 4. Анализ композиционного бета-разнообразия ---\n")

# CLR трансформация метаболома
fecal_clr_comp <- as.matrix(clr(fecal_num + 1))
aitchison_dist <- dist(fecal_clr_comp, method = "euclidean")

lipid_clr_comp <- as.matrix(clr(lipid_num + 1))
lipid_dist <- dist(lipid_clr_comp, method = "euclidean")

serum_clr_comp <- as.matrix(clr(serum_num + 1))
serum_dist <- dist(serum_clr_comp, method = "euclidean")

# PERMANOVA
set.seed(123)
permanova_res <- adonis2(aitchison_dist ~ y_model, permutations = 999)
permanova_res <- adonis2(lipid_dist ~ y_model, permutations = 999)
permanova_res <- adonis2(serum_dist ~ y_model, permutations = 999)
print(permanova_res)

p_val  <- permanova_res$`Pr(>F)`[1]
r2_val <- permanova_res$R2[1]

# PCoA
pcoa_res <- cmdscale(aitchison_dist, k = 2, eig = TRUE)
pcoa_res <- cmdscale(lipid_dist, k = 2, eig = TRUE)
pcoa_res <- cmdscale(serum_dist, k = 2, eig = TRUE)
variance_explained <- round(100 * pcoa_res$eig / sum(pcoa_res$eig), 1)

pcoa_df <- data.frame(PCoA1 = pcoa_res$points[, 1], PCoA2 = pcoa_res$points[, 2], Diagnosis = y_model)

# --- Шаг 4.1: Тест PERMDISP (Однородность дисперсий внутри групп) ---
cat("\n--- 4.1. Тест PERMDISP (Дисперсия внутри групп CD и UC) ---\n")

# Рассчитываем расстояния от центроидов групп
disp_res <- betadisper(aitchison_dist, group = y_model)
disp_res <- betadisper(lipid_dist, group = y_model)
disp_res <- betadisper(serum_dist, group = y_model)

# Статистический тест перестановками (999 перестановок)
set.seed(123)
permdisp_test <- permutest(disp_res, permutations = 999)
print(permdisp_test)

permdisp_pval <- permdisp_test$tab$`Pr(>F)`[1]
cat("PERMDISP p-value:", round(permdisp_pval, 4), "\n")

# Дополнительно: График дистанций до центроидов (boxplot)
boxplot(disp_res, main = "Dispersions from Centroids (PERMDISP)", 
        col = c("darkred", "darkblue"), xlab = "Diagnosis")


# --- Шаг 4.2: Попарный (Pairwise) PERMANOVA с поправкой ---
cat("\n--- 4.2. Попарный PERMANOVA (Сравнение групп) ---\n")

# Функция для попарного PERMANOVA (актуально для >2 групп, но применима и для двух)
pairwise_permanova <- function(dist_matrix, group_factor, padj_method = "BH") {
  co <- combn(unique(as.character(group_factor)), 2)
  pairs <- c()
  f_stats <- c()
  r2_scores <- c()
  p_values <- c()
  
  for(i in 1:ncol(co)) {
    # Выделяем маску для текущей пары групп
    mask <- group_factor %in% co[, i]
    sub_dist <- as.dist(as.matrix(dist_matrix)[mask, mask])
    sub_group <- factor(group_factor[mask])
    
    # Запускаем adonis2 локально
    set.seed(123)
    ad <- adonis2(sub_dist ~ sub_group, permutations = 999)
    
    pairs <- c(pairs, paste(co[1, i], "vs", co[2, i]))
    f_stats <- c(f_stats, ad$F[1])
    r2_scores <- c(r2_scores, ad$R2[1])
    p_values <- c(p_values, ad$`Pr(>F)`[1])
  }
  
  # Корректируем множественные сравнения по выбранному методу (по умолчанию Бенджамини-Хохберг)
  p_adjusted <- p.adjust(p_values, method = padj_method)
  
  results <- data.frame(
    Comparison = pairs,
    F_value = f_stats,
    R2 = r2_scores,
    p_value = p_values,
    p_adj = p_adjusted,
    Method = padj_method
  )
  return(results)
}

# Запуск попарного сравнения
pairwise_res <- pairwise_permanova(aitchison_dist, y_model, padj_method = "BH")
pairwise_res <- pairwise_permanova(lipid_dist, y_model, padj_method = "BH")
pairwise_res <- pairwise_permanova(serum_dist, y_model, padj_method = "BH")
print(pairwise_res)


# --- Шаг 4.3: Визуализация PCoA с нанесением всех статистических метрик ---
cat("\n--- 4.3. Визуализация PCoA с результатами тестов ---\n")

# Извлекаем финальные значения для вывода на график
p_permanova <- round(permanova_res$`Pr(>F)`[1], 4)
r2_permanova <- round(permanova_res$R2[1] * 100, 1)
p_permdisp  <- round(permdisp_pval, 4)
padj_pairwise <- round(pairwise_res$p_adj[1], 4) # Для пары CD vs UC

# Настройка текстовых строк для графика
stat_text <- c(
  paste0("Global PERMANOVA p = ", p_permanova, " (R² = ", r2_permanova, "%)"),
  paste0("PERMDISP (dispersion) p = ", p_permdisp),
  paste0("Pairwise p-adj (CD vs UC) = ", padj_pairwise)
)





pdf_output_path <- "Desktop/WORK/gut/1_stage/R/serum.pdf"
pdf(file = pdf_output_path, width = 8, height = 5.5, bg = "transparent")

par(
  bg = "transparent", 
  mar = c(4.5, 4.5, 4, 2), # Немного уменьшили внешние поля, так как текст теперь ближе
  mgp = c(2.2, 0.8, 0)           # 2.2 — расстояние названия, 0.8 — цифр, 0 — самой линии оси
) 

# --- Расчет динамических лимитов с запасом 25% для вмещения эллипсов ---
xlim_range <- range(pcoa_df$PCoA1)
ylim_range <- range(pcoa_df$PCoA2)

# Добавляем отступы во все стороны
xlim_padded <- xlim_range + c(-0.25, 0.25) * diff(xlim_range)
ylim_padded <- ylim_range + c(-0.25, 0.25) * diff(ylim_range)

# Отрисовка графика
plot(pcoa_df$PCoA1, pcoa_df$PCoA2, 
     col = ifelse(pcoa_df$Diagnosis == "Class_CD", "#D62828", "#7209B7"), 
     pch = 19, cex = 1.5,
     xlim = xlim_padded, 
     ylim = ylim_padded,
     cex.lab = 1.4,      # Увеличенный размер названий осей
     cex.axis = 1.2,     # Увеличенный размер цифр на осях
     cex.main = 1.4,     # Увеличенный размер заголовка
     xlab = paste0("PCoA 1 (", variance_explained[1], "%)"), 
     ylab = paste0("PCoA 2 (", variance_explained[2], "%)"),
     main = "PCoA Serum Metabolome (Aitchison distance)")

# Добавление эллипсов
ordiellipse(pcoa_res$points, y_model, kind = "ehull", 
            col = c("#D62828", "#7209B7"), lwd = 2, lty = 2)

# Добавление легенды групп
legend("topright", legend = c("CD", "UC"), 
       col = c("#D62828", "#7209B7"), pch = 19, bty = "n", 
       cex = 1.2) 


# Добавление блока со статистикой (в левый нижний угол)
legend("bottomright", legend = stat_text, 
       bty = "o", box.col = "gray", bg = "white", cex = 0.65, title = "Statistics:")

dev.off()














################histology



# ============================================================================
# АНАЛИЗ: FECAL МЕТАБОЛОМ КАК ПРЕДИКТОР ГИСТОЛОГИЧЕСКИХ ПРИЗНАКОВ
# (ДЛЯ ГОТОВЫХ ОТФИЛЬТРОВАННЫХ ФАЙЛОВ)
# ============================================================================

# --- Шаг 1: Загрузка данных ---
fecal_data <- read.table("Desktop/WORK/gut/1_stage/R/fecal_for_histology.txt", 
                         header = TRUE, row.names = 1, sep = "\t", 
                         check.names = FALSE)

histology_data <- read.table("Desktop/WORK/gut/1_stage/R/histology.txt", 
                             header = TRUE, row.names = 1, sep = "\t", 
                             check.names = FALSE)

cat("\n=========================================")
cat("\nFECAL МЕТАБОЛОМ КАК ПРЕДИКТОР ГИСТОЛОГИИ")
cat("\n=========================================\n")

cat("\n--- Размеры данных ---\n")
cat("  Фекальные метаболиты:", dim(fecal_data), "\n")
cat("  Гистология:", dim(histology_data), "\n")

# --- Шаг 2: Предобработка фекальных данных ---

cat("\n--- Предобработка фекальных метаболитов ---\n")

# Проверяем, что данные уже отфильтрованы
if(any(is.na(fecal_data))) {
  cat("  ⚠️ Обнаружены NA, заменяем на 0\n")
  fecal_data[is.na(fecal_data)] <- 0
}

# Логарифмическая трансформация и масштабирование
fecal_log <- log2(fecal_data + 1)
fecal_scaled <- scale(fecal_log)

# Удаляем признаки с нулевой дисперсией
variances <- apply(fecal_scaled, 2, var, na.rm = TRUE)
fecal_ready <- fecal_scaled[, variances > 0 & !is.na(variances), drop = FALSE]

cat("  Размер после предобработки:", dim(fecal_ready), "\n")

# --- Шаг 3: Подготовка гистологических данных ---

cat("\n--- Подготовка гистологических данных ---\n")

# Проверяем, что данные бинарные (0/1)
histology_matrix <- as.matrix(histology_data)

# Убеждаемся, что все значения 0 или 1
histology_matrix[histology_matrix != 0 & histology_matrix != 1] <- 0

cat("  Размер матрицы гистологии:", dim(histology_matrix), "\n")
cat("  Количество признаков:", ncol(histology_matrix), "\n")

# Статистика по каждому признаку
for(j in 1:ncol(histology_matrix)) {
  n_pos <- sum(histology_matrix[, j] == 1, na.rm = TRUE)
  n_neg <- sum(histology_matrix[, j] == 0, na.rm = TRUE)
  cat("    ", colnames(histology_matrix)[j], ":", n_pos, "положительных,", n_neg, "отрицательных\n")
}

# --- Шаг 4: Функция для проверки предсказательной способности ---

test_predictor <- function(fecal_data, y_binary, feature_name) {
  
  # Создаем пустой результат
  result <- data.frame(
    Feature = feature_name,
    N_Positive = NA_integer_,
    N_Negative = NA_integer_,
    Wilcoxon_p = NA_real_,
    RF_Accuracy = NA_real_,
    PLSDA_Accuracy = NA_real_,
    stringsAsFactors = FALSE
  )
  
  # Проверяем y_binary
  if(is.null(y_binary) || length(y_binary) == 0) {
    return(result)
  }
  
  # Убираем NA
  valid_idx <- !is.na(y_binary)
  if(sum(valid_idx) < 3) {
    return(result)
  }
  
  y_clean <- y_binary[valid_idx]
  X_clean <- fecal_data[valid_idx, , drop = FALSE]
  
  # Проверяем, что есть оба класса
  unique_vals <- unique(y_clean)
  if(length(unique_vals) < 2) {
    return(result)
  }
  
  # Проверяем X_clean
  if(ncol(X_clean) == 0) {
    return(result)
  }
  
  # Заполняем N_Positive и N_Negative
  result$N_Positive <- sum(y_clean == 1, na.rm = TRUE)
  result$N_Negative <- sum(y_clean == 0, na.rm = TRUE)
  
  # --- Wilcoxon ---
  tryCatch({
    p_values <- apply(X_clean, 2, function(x) {
      if(var(x, na.rm = TRUE) == 0 || length(unique(x)) < 2) return(NA)
      test <- wilcox.test(x ~ y_clean)
      return(test$p.value)
    })
    
    min_p <- min(p_values, na.rm = TRUE)
    if(!is.infinite(min_p) && !is.na(min_p)) {
      adj_p <- min_p * ncol(X_clean)
      result$Wilcoxon_p <- min(adj_p, 1)
    }
  }, error = function(e) {
    # Оставляем NA
  })
  
  # --- Random Forest ---
  tryCatch({
    library(randomForest)
    set.seed(123)
    
    if(ncol(X_clean) > 50) {
      vars <- apply(X_clean, 2, var, na.rm = TRUE)
      top_vars <- order(vars, decreasing = TRUE)[1:min(50, length(vars))]
      X_clean_rf <- X_clean[, top_vars, drop = FALSE]
    } else {
      X_clean_rf <- X_clean
    }
    
    if(ncol(X_clean_rf) > 1) {
      rf_model <- randomForest(x = X_clean_rf, y = as.factor(y_clean), 
                               ntree = 100, importance = TRUE)
      result$RF_Accuracy <- 1 - rf_model$err.rate[nrow(rf_model$err.rate), "OOB"]
    }
  }, error = function(e) {
    # Оставляем NA
  })
  
  # --- PLS-DA ---
  tryCatch({
    library(mixOmics)
    set.seed(123)
    
    if(ncol(X_clean) > 50) {
      vars <- apply(X_clean, 2, var, na.rm = TRUE)
      top_vars <- order(vars, decreasing = TRUE)[1:min(50, length(vars))]
      X_clean_pls <- X_clean[, top_vars, drop = FALSE]
    } else {
      X_clean_pls <- X_clean
    }
    
    if(ncol(X_clean_pls) > 1) {
      plsda_model <- plsda(X_clean_pls, as.factor(y_clean), ncomp = 2)
      perf_plsda <- perf(plsda_model, validation = 'Mfold', folds = 5, nrepeat = 3)
      error_rate <- perf_plsda$error.rate$max.dist[1]
      result$PLSDA_Accuracy <- 1 - error_rate
    }
  }, error = function(e) {
    # Оставляем NA
  })
  
  return(result)
}

# --- Шаг 5: Проверка предсказательной способности для каждого признака ---

cat("\n--- Проверка предсказательной способности ---\n")

results_all <- data.frame()

for(j in 1:ncol(histology_matrix)) {
  feature_name <- colnames(histology_matrix)[j]
  
  cat("\n  Обработка:", feature_name, "\n")
  
  y_binary <- histology_matrix[, j]
  n_positive <- sum(y_binary == 1, na.rm = TRUE)
  n_negative <- sum(y_binary == 0, na.rm = TRUE)
  
  cat("    Положительных:", n_positive, ", Отрицательных:", n_negative, "\n")
  
  if(n_positive < 3 || n_negative < 3) {
    cat("    ⚠️ Слишком мало данных, пропускаем\n")
    next
  }
  
  result <- test_predictor(fecal_ready, y_binary, feature_name)
  results_all <- rbind(results_all, result)
}







# --- Шаг 6: Результаты (максимально просто) ---

cat("\n=========================================")
cat("\nРЕЗУЛЬТАТЫ")
cat("\n=========================================\n")

# Сортируем по p-value
results_all <- results_all[order(results_all$Wilcoxon_p), ]

# ---- ВЫВОД ТАБЛИЦЫ ----
cat("\n📊 Сводная таблица:\n")
cat("------------------------------------------------------------\n")
cat(sprintf("%-40s %8s %8s %12s %12s\n", 
            "Признак", "N_pos", "N_neg", "p-value", "RF_Acc"))
cat("------------------------------------------------------------\n")

for(i in 1:nrow(results_all)) {
  cat(sprintf("%-40s %8d %8d %12.4f %12.4f\n",
              substr(results_all$Feature[i], 1, 40),
              results_all$N_Positive[i],
              results_all$N_Negative[i],
              results_all$Wilcoxon_p[i],
              results_all$RF_Accuracy[i]))
}
cat("------------------------------------------------------------\n")

# ---- БЛИЖАЙШИЕ К ЗНАЧИМОСТИ ----
cat("\n🔍 Ближайшие к значимости (p < 0.1):\n")
near_sig <- results_all[results_all$Wilcoxon_p < 0.1 & !is.na(results_all$Wilcoxon_p), ]
for(i in 1:nrow(near_sig)) {
  cat(sprintf("  %s: p = %.4f, RF_Acc = %.3f\n",
              near_sig$Feature[i],
              near_sig$Wilcoxon_p[i],
              near_sig$RF_Accuracy[i]))
}
if(nrow(near_sig) == 0) cat("  Нет признаков с p < 0.1\n")

# ---- ЛУЧШИЕ ПО ТОЧНОСТИ ----
cat("\n🏆 Лучшие по точности Random Forest:\n")
best_rf <- results_all[!is.na(results_all$RF_Accuracy), ]
best_rf <- best_rf[order(-best_rf$RF_Accuracy), ]
best_rf <- best_rf[1:min(3, nrow(best_rf)), ]
for(i in 1:nrow(best_rf)) {
  cat(sprintf("  %s: RF_Acc = %.3f, p = %.4f\n",
              best_rf$Feature[i],
              best_rf$RF_Accuracy[i],
              best_rf$Wilcoxon_p[i]))
}
if(nrow(best_rf) == 0) cat("  Нет данных\n")

# ---- ЗНАЧИМЫЕ РЕЗУЛЬТАТЫ ----
cat("\n✅ Значимые результаты (p < 0.05):\n")
sig <- results_all[results_all$Wilcoxon_p < 0.05 & !is.na(results_all$Wilcoxon_p), ]
for(i in 1:nrow(sig)) {
  cat(sprintf("  %s: p = %.4f\n", sig$Feature[i], sig$Wilcoxon_p[i]))
}
if(nrow(sig) == 0) cat("  ❌ Нет значимых результатов\n")

# ---- ИТОГ ----
cat("\n📌 ИТОГ:\n")
cat(sprintf("  Всего признаков: %d\n", nrow(results_all)))
cat(sprintf("  Значимых (p < 0.05): %d\n", sum(results_all$Wilcoxon_p < 0.05, na.rm = TRUE)))
cat(sprintf("  На границе (p < 0.1): %d\n", sum(results_all$Wilcoxon_p < 0.1 & results_all$Wilcoxon_p >= 0.05, na.rm = TRUE)))
cat(sprintf("  Средняя точность RF: %.3f\n", mean(results_all$RF_Accuracy, na.rm = TRUE)))

# Сохраняем
write.table(results_all, "histology_predictors_results.tsv", 
            sep = "\t", row.names = FALSE, quote = FALSE)
cat("\n✅ Результаты сохранены в: histology_predictors_results.tsv\n")









# --- Шаг 6: Результаты ---

# ============================================================================
# ВИЗУАЛИЗАЦИЯ: РИСУНКИ ПО ОДНОМУ
# ============================================================================

library(ggplot2)

cat("\n=========================================")
cat("\nСОЗДАНИЕ РИСУНКОВ")
cat("\n=========================================\n")

# --- Подготовка данных для графиков ---
plot_data <- results_all
plot_data$log_p <- -log10(plot_data$Wilcoxon_p)
plot_data$Significant <- plot_data$Wilcoxon_p < 0.05
plot_data$Near_Significant <- plot_data$Wilcoxon_p < 0.06 & plot_data$Wilcoxon_p >= 0.05

plot_data$Group <- ifelse(plot_data$Significant, "Significant (p < 0.05)",
                          ifelse(plot_data$Near_Significant, "Near Significant (p < 0.06)",
                                 "Not Significant"))

# Укорачиваем названия
plot_data$Feature_Short <- gsub("terminal_ileum_", "TI_", plot_data$Feature)
plot_data$Feature_Short <- gsub("cecal_pole_", "CP_", plot_data$Feature_Short)
plot_data$Feature_Short <- gsub("_", " ", plot_data$Feature_Short)


  plot_data_combined <- plot_data
  plot_data_combined$Accuracy_Category <- ifelse(plot_data_combined$RF_Accuracy > 0.7, "High (>0.7)", 
                                                 ifelse(plot_data_combined$RF_Accuracy > 0.5, "Medium (>0.5)", "Low"))
  

  
  plot_data_combined$color_group <- ifelse(grepl("CP 3", plot_data_combined$Feature_Short, ignore.case = TRUE), 
                                           "red", "gray")
  
  pdf("Desktop/WORK/gut/1_stage/R/histology_combined_plot.pdf", width = 8, height = 5.5)
  
  ggplot(plot_data_combined, aes(x = Wilcoxon_p, 
                                 y = RF_Accuracy,
                                 size = N_Positive + N_Negative,
                                 color = color_group
  )) +
    geom_point(alpha = 0.8) +
    geom_vline(xintercept = 0.05, linetype = "dashed", color = "red", size = 0.8) +
    geom_hline(yintercept = 0.85, linetype = "dashed", color = "red", size = 0.8) +
    labs(title = "Fecal metabolome VS. Histology",
         subtitle = "",
         x = "P-value",
         y = "Accuracy (Random Forest)") +
    theme_minimal() +
    theme(
          axis.text = element_text(size = 17),
          plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
          axis.text.y = element_text(size = 14),
          axis.text.x = element_text(size = 14),
          legend.position = "none") +
    xlim(0, 1) +
    scale_color_identity()  # используем цвета как есть
  
  dev.off()
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  ###################################Drug
  
  
  # ============================================================================
  # АНАЛИЗ: FECAL МЕТАБОЛОМ КАК ПРЕДИКТОР ТЕРАПИИ
  # ============================================================================
  
  # --- Шаг 1: Загрузка данных ---
  fecal_data <- read.table("Desktop/WORK/gut/1_stage/R/fecal_for_drug.txt", 
                           header = TRUE, row.names = 1, sep = "\t", 
                           check.names = FALSE)
  
  drug_data <- read.table("Desktop/WORK/gut/1_stage/R/drug.txt", 
                          header = TRUE, row.names = 1, sep = "\t", 
                          check.names = FALSE)
  
  cat("\n=========================================")
  cat("\nFECAL МЕТАБОЛОМ КАК ПРЕДИКТОР ТЕРАПИИ")
  cat("\n=========================================\n")
  
  cat("\n--- Размеры данных ---\n")
  cat("  Фекальные метаболиты:", dim(fecal_data), "\n")
  cat("  Лекарства:", dim(drug_data), "\n")
  
  # --- Шаг 2: Предобработка фекальных данных ---
  
  cat("\n--- Предобработка фекальных метаболитов ---\n")
  
  # Логарифмическая трансформация и масштабирование
  fecal_log <- log2(fecal_data + 1)
  fecal_scaled <- scale(fecal_log)
  
  # Удаляем признаки с нулевой дисперсией
  variances <- apply(fecal_scaled, 2, var, na.rm = TRUE)
  fecal_ready <- fecal_scaled[, variances > 0 & !is.na(variances), drop = FALSE]
  
  cat("  Размер после предобработки:", dim(fecal_ready), "\n")
  
  # --- Шаг 3: Подготовка данных о лекарствах ---
  
  cat("\n--- Подготовка данных о лекарствах ---\n")
  
  drug_matrix <- as.matrix(drug_data)
  drug_matrix[drug_matrix != 0 & drug_matrix != 1] <- 0
  
  cat("  Размер матрицы лекарств:", dim(drug_matrix), "\n")
  cat("  Количество препаратов:", ncol(drug_matrix), "\n")
  
  # Статистика по каждому препарату
  cat("\n  Статистика по препаратам:\n")
  for(j in 1:ncol(drug_matrix)) {
    n_pos <- sum(drug_matrix[, j] == 1, na.rm = TRUE)
    n_neg <- sum(drug_matrix[, j] == 0, na.rm = TRUE)
    cat("    ", colnames(drug_matrix)[j], ":", n_pos, "принимают,", n_neg, "не принимают\n")
  }
  
  # --- Шаг 4: Функция для проверки предсказательной способности ---
  
  test_predictor <- function(fecal_data, y_binary, feature_name) {
    
    result <- data.frame(
      Feature = feature_name,
      N_Positive = NA_integer_,
      N_Negative = NA_integer_,
      Wilcoxon_p = NA_real_,
      RF_Accuracy = NA_real_,
      PLSDA_Accuracy = NA_real_,
      stringsAsFactors = FALSE
    )
    
    if(is.null(y_binary) || length(y_binary) == 0) return(result)
    
    valid_idx <- !is.na(y_binary)
    if(sum(valid_idx) < 3) return(result)
    
    y_clean <- y_binary[valid_idx]
    X_clean <- fecal_data[valid_idx, , drop = FALSE]
    
    unique_vals <- unique(y_clean)
    if(length(unique_vals) < 2) return(result)
    if(ncol(X_clean) == 0) return(result)
    
    result$N_Positive <- sum(y_clean == 1, na.rm = TRUE)
    result$N_Negative <- sum(y_clean == 0, na.rm = TRUE)
    
    # --- Wilcoxon ---
    tryCatch({
      p_values <- apply(X_clean, 2, function(x) {
        if(var(x, na.rm = TRUE) == 0 || length(unique(x)) < 2) return(NA)
        test <- wilcox.test(x ~ y_clean)
        return(test$p.value)
      })
      
      min_p <- min(p_values, na.rm = TRUE)
      if(!is.infinite(min_p) && !is.na(min_p)) {
        adj_p <- min_p * ncol(X_clean)
        result$Wilcoxon_p <- min(adj_p, 1)
      }
    }, error = function(e) {})
    
    # --- Random Forest ---
    tryCatch({
      library(randomForest)
      set.seed(123)
      
      if(ncol(X_clean) > 50) {
        vars <- apply(X_clean, 2, var, na.rm = TRUE)
        top_vars <- order(vars, decreasing = TRUE)[1:min(50, length(vars))]
        X_clean_rf <- X_clean[, top_vars, drop = FALSE]
      } else {
        X_clean_rf <- X_clean
      }
      
      if(ncol(X_clean_rf) > 1) {
        rf_model <- randomForest(x = X_clean_rf, y = as.factor(y_clean), 
                                 ntree = 100, importance = TRUE)
        result$RF_Accuracy <- 1 - rf_model$err.rate[nrow(rf_model$err.rate), "OOB"]
      }
    }, error = function(e) {})
    
    # --- PLS-DA ---
    tryCatch({
      library(mixOmics)
      set.seed(123)
      
      if(ncol(X_clean) > 50) {
        vars <- apply(X_clean, 2, var, na.rm = TRUE)
        top_vars <- order(vars, decreasing = TRUE)[1:min(50, length(vars))]
        X_clean_pls <- X_clean[, top_vars, drop = FALSE]
      } else {
        X_clean_pls <- X_clean
      }
      
      if(ncol(X_clean_pls) > 1) {
        plsda_model <- plsda(X_clean_pls, as.factor(y_clean), ncomp = 2)
        perf_plsda <- perf(plsda_model, validation = 'Mfold', folds = 5, nrepeat = 3)
        error_rate <- perf_plsda$error.rate$max.dist[1]
        result$PLSDA_Accuracy <- 1 - error_rate
      }
    }, error = function(e) {})
    
    return(result)
  }
  
  # --- Шаг 5: Проверка предсказательной способности для каждого препарата ---
  
  cat("\n--- Проверка предсказательной способности ---\n")
  
  results_all <- data.frame()
  
  for(j in 1:ncol(drug_matrix)) {
    feature_name <- colnames(drug_matrix)[j]
    
    cat("\n  Обработка:", feature_name, "\n")
    
    y_binary <- drug_matrix[, j]
    n_positive <- sum(y_binary == 1, na.rm = TRUE)
    n_negative <- sum(y_binary == 0, na.rm = TRUE)
    
    cat("    Принимают:", n_positive, ", Не принимают:", n_negative, "\n")
    
    if(n_positive < 3 || n_negative < 3) {
      cat("    ⚠️ Слишком мало данных, пропускаем\n")
      next
    }
    
    result <- test_predictor(fecal_ready, y_binary, feature_name)
    results_all <- rbind(results_all, result)
  }
  
  # --- Шаг 6: Результаты ---
  
  cat("\n=========================================")
  cat("\nРЕЗУЛЬТАТЫ")
  cat("\n=========================================\n")
  
  if(nrow(results_all) > 0) {
    
    results_all <- results_all[order(results_all$Wilcoxon_p), ]
    
    cat("\n📊 Сводная таблица:\n")
    cat("------------------------------------------------------------\n")
    cat(sprintf("%-25s %8s %8s %12s %12s\n", 
                "Препарат", "Принимают", "Не принимают", "p-value", "RF_Acc"))
    cat("------------------------------------------------------------\n")
    
    for(i in 1:nrow(results_all)) {
      cat(sprintf("%-25s %8d %8d %12.4f %12.4f\n",
                  substr(results_all$Feature[i], 1, 25),
                  results_all$N_Positive[i],
                  results_all$N_Negative[i],
                  results_all$Wilcoxon_p[i],
                  results_all$RF_Accuracy[i]))
    }
    cat("------------------------------------------------------------\n")
    
    # ---- БЛИЖАЙШИЕ К ЗНАЧИМОСТИ ----
    cat("\n🔍 Ближайшие к значимости (p < 0.1):\n")
    near_sig <- results_all[results_all$Wilcoxon_p < 0.1 & !is.na(results_all$Wilcoxon_p), ]
    if(nrow(near_sig) > 0) {
      for(i in 1:nrow(near_sig)) {
        cat(sprintf("  %s: p = %.4f, RF_Acc = %.3f\n",
                    near_sig$Feature[i],
                    near_sig$Wilcoxon_p[i],
                    near_sig$RF_Accuracy[i]))
      }
    } else {
      cat("  Нет препаратов с p < 0.1\n")
    }
    
    # ---- ЛУЧШИЕ ПО ТОЧНОСТИ ----
    cat("\n🏆 Лучшие по точности Random Forest:\n")
    best_rf <- results_all[!is.na(results_all$RF_Accuracy), ]
    best_rf <- best_rf[order(-best_rf$RF_Accuracy), ]
    best_rf <- best_rf[1:min(3, nrow(best_rf)), ]
    for(i in 1:nrow(best_rf)) {
      cat(sprintf("  %s: RF_Acc = %.3f, p = %.4f\n",
                  best_rf$Feature[i],
                  best_rf$RF_Accuracy[i],
                  best_rf$Wilcoxon_p[i]))
    }
    
    # ---- ЗНАЧИМЫЕ РЕЗУЛЬТАТЫ ----
    cat("\n✅ Значимые результаты (p < 0.05):\n")
    sig <- results_all[results_all$Wilcoxon_p < 0.05 & !is.na(results_all$Wilcoxon_p), ]
    if(nrow(sig) > 0) {
      for(i in 1:nrow(sig)) {
        cat(sprintf("  %s: p = %.4f\n", sig$Feature[i], sig$Wilcoxon_p[i]))
      }
    } else {
      cat("  ❌ Нет значимых результатов\n")
    }
    
    # ---- ИТОГ ----
    cat("\n📌 ИТОГ:\n")
    cat(sprintf("  Всего препаратов: %d\n", nrow(results_all)))
    cat(sprintf("  Значимых (p < 0.05): %d\n", sum(results_all$Wilcoxon_p < 0.05, na.rm = TRUE)))
    cat(sprintf("  На границе (p < 0.1): %d\n", sum(results_all$Wilcoxon_p < 0.1 & results_all$Wilcoxon_p >= 0.05, na.rm = TRUE)))
    cat(sprintf("  Средняя точность RF: %.3f\n", mean(results_all$RF_Accuracy, na.rm = TRUE)))
    
    # Сохраняем
    write.table(results_all, "drug_predictors_results.tsv", 
                sep = "\t", row.names = FALSE, quote = FALSE)
    cat("\n✅ Результаты сохранены в: drug_predictors_results.tsv\n")
    
  } else {
    cat("\n❌ Нет результатов для анализа\n")
  }
  
  cat("\n=========================================")
  cat("\nАНАЛИЗ ЗАВЕРШЕН!")
  cat("\n=========================================\n")
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  # ============================================================================
  # АНАЛИЗ: FECAL МЕТАБОЛОМ КАК ПРЕДИКТОР КОМБИНАТОРНОЙ ТЕРАПИИ
  # ============================================================================
  
  # --- Шаг 1: Загрузка данных ---
  fecal_data <- read.table("Desktop/WORK/gut/1_stage/R/fecal_for_drug.txt", 
                           header = TRUE, row.names = 1, sep = "\t", 
                           check.names = FALSE)
  
  drug_data <- read.table("Desktop/WORK/gut/1_stage/R/drug.txt", 
                          header = TRUE, row.names = 1, sep = "\t", 
                          check.names = FALSE)
  
  cat("\n=========================================")
  cat("\nFECAL МЕТАБОЛОМ КАК ПРЕДИКТОР КОМБИНАТОРНОЙ ТЕРАПИИ")
  cat("\n=========================================\n")
  
  # --- Шаг 2: Предобработка ---
  
  cat("\n--- Предобработка фекальных метаболитов ---\n")
  
  fecal_log <- log2(fecal_data + 1)
  fecal_scaled <- scale(fecal_log)
  
  variances <- apply(fecal_scaled, 2, var, na.rm = TRUE)
  fecal_ready <- fecal_scaled[, variances > 0 & !is.na(variances), drop = FALSE]
  
  cat("  Размер фекальных данных:", dim(fecal_ready), "\n")
  
  # --- Шаг 3: Создание комбинаторных признаков ---
  
  cat("\n--- Создание комбинаций препаратов ---\n")
  
  drug_matrix <- as.matrix(drug_data)
  drug_matrix[drug_matrix != 0 & drug_matrix != 1] <- 0
  
  drug_names <- colnames(drug_matrix)
  n_drugs <- ncol(drug_matrix)
  
  cat("  Препараты:", paste(drug_names, collapse = ", "), "\n")
  
  # 3.1 Отдельные препараты
  cat("\n  Отдельные препараты:\n")
  for(j in 1:n_drugs) {
    n_pos <- sum(drug_matrix[, j] == 1, na.rm = TRUE)
    cat("    ", drug_names[j], ":", n_pos, "пациентов\n")
  }
  
  # 3.2 Парные комбинации
  cat("\n  Парные комбинации:\n")
  pair_combinations <- combn(drug_names, 2, simplify = FALSE)
  pair_names <- sapply(pair_combinations, function(x) paste(x, collapse = "+"))
  
  for(k in 1:length(pair_combinations)) {
    pair <- pair_combinations[[k]]
    combo <- rowSums(drug_matrix[, pair]) == 2
    n_combo <- sum(combo, na.rm = TRUE)
    if(n_combo >= 3) {
      cat("    ", pair_names[k], ":", n_combo, "пациентов\n")
    }
  }
  
  # 3.3 Тройные комбинации
  cat("\n  Тройные комбинации:\n")
  triple_combinations <- combn(drug_names, 3, simplify = FALSE)
  triple_names <- sapply(triple_combinations, function(x) paste(x, collapse = "+"))
  
  for(k in 1:length(triple_combinations)) {
    triple <- triple_combinations[[k]]
    combo <- rowSums(drug_matrix[, triple]) == 3
    n_combo <- sum(combo, na.rm = TRUE)
    if(n_combo >= 3) {
      cat("    ", triple_names[k], ":", n_combo, "пациентов\n")
    }
  }
  
  # --- Шаг 4: Функция для проверки предсказательной способности ---
  
  test_predictor <- function(fecal_data, y_binary, feature_name) {
    
    result <- data.frame(
      Feature = feature_name,
      N_Positive = NA_integer_,
      N_Negative = NA_integer_,
      Wilcoxon_p = NA_real_,
      RF_Accuracy = NA_real_,
      PLSDA_Accuracy = NA_real_,
      stringsAsFactors = FALSE
    )
    
    if(is.null(y_binary) || length(y_binary) == 0) return(result)
    
    valid_idx <- !is.na(y_binary)
    if(sum(valid_idx) < 3) return(result)
    
    y_clean <- y_binary[valid_idx]
    X_clean <- fecal_data[valid_idx, , drop = FALSE]
    
    unique_vals <- unique(y_clean)
    if(length(unique_vals) < 2) return(result)
    if(ncol(X_clean) == 0) return(result)
    
    result$N_Positive <- sum(y_clean == 1, na.rm = TRUE)
    result$N_Negative <- sum(y_clean == 0, na.rm = TRUE)
    
    # --- Wilcoxon ---
    tryCatch({
      p_values <- apply(X_clean, 2, function(x) {
        if(var(x, na.rm = TRUE) == 0 || length(unique(x)) < 2) return(NA)
        test <- wilcox.test(x ~ y_clean)
        return(test$p.value)
      })
      
      min_p <- min(p_values, na.rm = TRUE)
      if(!is.infinite(min_p) && !is.na(min_p)) {
        adj_p <- min_p * ncol(X_clean)
        result$Wilcoxon_p <- min(adj_p, 1)
      }
    }, error = function(e) {})
    
    # --- Random Forest ---
    tryCatch({
      library(randomForest)
      set.seed(123)
      
      if(ncol(X_clean) > 50) {
        vars <- apply(X_clean, 2, var, na.rm = TRUE)
        top_vars <- order(vars, decreasing = TRUE)[1:min(50, length(vars))]
        X_clean_rf <- X_clean[, top_vars, drop = FALSE]
      } else {
        X_clean_rf <- X_clean
      }
      
      if(ncol(X_clean_rf) > 1) {
        rf_model <- randomForest(x = X_clean_rf, y = as.factor(y_clean), 
                                 ntree = 100, importance = TRUE)
        result$RF_Accuracy <- 1 - rf_model$err.rate[nrow(rf_model$err.rate), "OOB"]
      }
    }, error = function(e) {})
    
    # --- PLS-DA ---
    tryCatch({
      library(mixOmics)
      set.seed(123)
      
      if(ncol(X_clean) > 50) {
        vars <- apply(X_clean, 2, var, na.rm = TRUE)
        top_vars <- order(vars, decreasing = TRUE)[1:min(50, length(vars))]
        X_clean_pls <- X_clean[, top_vars, drop = FALSE]
      } else {
        X_clean_pls <- X_clean
      }
      
      if(ncol(X_clean_pls) > 1) {
        plsda_model <- plsda(X_clean_pls, as.factor(y_clean), ncomp = 2)
        perf_plsda <- perf(plsda_model, validation = 'Mfold', folds = 5, nrepeat = 3)
        error_rate <- perf_plsda$error.rate$max.dist[1]
        result$PLSDA_Accuracy <- 1 - error_rate
      }
    }, error = function(e) {})
    
    return(result)
  }
  
  # --- Шаг 5: Сбор всех комбинаций для анализа ---
  
  cat("\n--- Сбор комбинаций для анализа ---\n")
  
  all_combinations <- list()
  
  # 5.1 Отдельные препараты
  for(j in 1:n_drugs) {
    name <- drug_names[j]
    all_combinations[[name]] <- list(
      name = name,
      y = drug_matrix[, j],
      type = "single"
    )
  }
  
  # 5.2 Парные комбинации
  for(k in 1:length(pair_combinations)) {
    pair <- pair_combinations[[k]]
    combo <- rowSums(drug_matrix[, pair]) == 2
    name <- paste(pair, collapse = "+")
    if(sum(combo, na.rm = TRUE) >= 3) {
      all_combinations[[name]] <- list(
        name = name,
        y = as.numeric(combo),
        type = "pair"
      )
    }
  }
  
  # 5.3 Тройные комбинации
  for(k in 1:length(triple_combinations)) {
    triple <- triple_combinations[[k]]
    combo <- rowSums(drug_matrix[, triple]) == 3
    name <- paste(triple, collapse = "+")
    if(sum(combo, na.rm = TRUE) >= 3) {
      all_combinations[[name]] <- list(
        name = name,
        y = as.numeric(combo),
        type = "triple"
      )
    }
  }
  
  cat("  Всего комбинаций для анализа:", length(all_combinations), "\n")
  
  # --- Шаг 6: Проверка предсказательной способности ---
  
  cat("\n--- Проверка предсказательной способности ---\n")
  
  results_all <- data.frame()
  
  for(combo_name in names(all_combinations)) {
    combo_info <- all_combinations[[combo_name]]
    
    cat("\n  Обработка:", combo_name, "(", combo_info$type, ")\n")
    
    y_binary <- combo_info$y
    n_positive <- sum(y_binary == 1, na.rm = TRUE)
    n_negative <- sum(y_binary == 0, na.rm = TRUE)
    
    cat("    Положительных:", n_positive, ", Отрицательных:", n_negative, "\n")
    
    if(n_positive < 3 || n_negative < 3) {
      cat("    ⚠️ Слишком мало данных, пропускаем\n")
      next
    }
    
    result <- test_predictor(fecal_ready, y_binary, combo_name)
    result$Type <- combo_info$type
    results_all <- rbind(results_all, result)
  }
  
  # --- Шаг 7: Результаты ---
  
  cat("\n=========================================")
  cat("\nРЕЗУЛЬТАТЫ")
  cat("\n=========================================\n")
  
  if(nrow(results_all) > 0) {
    
    results_all <- results_all[order(results_all$Wilcoxon_p), ]
    
    cat("\n📊 Сводная таблица:\n")
    cat("--------------------------------------------------------------------\n")
    cat(sprintf("%-30s %8s %8s %12s %12s %10s\n", 
                "Комбинация", "Принимают", "Не принимают", "p-value", "RF_Acc", "Тип"))
    cat("--------------------------------------------------------------------\n")
    
    for(i in 1:nrow(results_all)) {
      cat(sprintf("%-30s %8d %8d %12.4f %12.4f %10s\n",
                  substr(results_all$Feature[i], 1, 30),
                  results_all$N_Positive[i],
                  results_all$N_Negative[i],
                  results_all$Wilcoxon_p[i],
                  results_all$RF_Accuracy[i],
                  results_all$Type[i]))
    }
    cat("--------------------------------------------------------------------\n")
    
    # ---- ЛУЧШИЕ ПО ТОЧНОСТИ ----
    cat("\n🏆 Лучшие по точности Random Forest:\n")
    best_rf <- results_all[!is.na(results_all$RF_Accuracy), ]
    best_rf <- best_rf[order(-best_rf$RF_Accuracy), ]
    best_rf <- best_rf[1:min(5, nrow(best_rf)), ]
    for(i in 1:nrow(best_rf)) {
      cat(sprintf("  %s: RF_Acc = %.3f, p = %.4f (%s)\n",
                  best_rf$Feature[i],
                  best_rf$RF_Accuracy[i],
                  best_rf$Wilcoxon_p[i],
                  best_rf$Type[i]))
    }
    
    # ---- БЛИЖАЙШИЕ К ЗНАЧИМОСТИ ----
    cat("\n🔍 Ближайшие к значимости (p < 0.1):\n")
    near_sig <- results_all[results_all$Wilcoxon_p < 0.1 & !is.na(results_all$Wilcoxon_p), ]
    if(nrow(near_sig) > 0) {
      for(i in 1:nrow(near_sig)) {
        cat(sprintf("  %s: p = %.4f, RF_Acc = %.3f (%s)\n",
                    near_sig$Feature[i],
                    near_sig$Wilcoxon_p[i],
                    near_sig$RF_Accuracy[i],
                    near_sig$Type[i]))
      }
    } else {
      cat("  Нет комбинаций с p < 0.1\n")
    }
    
    # ---- ЗНАЧИМЫЕ РЕЗУЛЬТАТЫ ----
    cat("\n✅ Значимые результаты (p < 0.05):\n")
    sig <- results_all[results_all$Wilcoxon_p < 0.05 & !is.na(results_all$Wilcoxon_p), ]
    if(nrow(sig) > 0) {
      for(i in 1:nrow(sig)) {
        cat(sprintf("  %s: p = %.4f (%s)\n", 
                    sig$Feature[i], 
                    sig$Wilcoxon_p[i],
                    sig$Type[i]))
      }
    } else {
      cat("  ❌ Нет значимых результатов\n")
    }
    
    # ---- ИТОГ ----
    cat("\n📌 ИТОГ:\n")
    cat(sprintf("  Всего комбинаций: %d\n", nrow(results_all)))
    cat(sprintf("  Значимых (p < 0.05): %d\n", sum(results_all$Wilcoxon_p < 0.05, na.rm = TRUE)))
    cat(sprintf("  На границе (p < 0.1): %d\n", sum(results_all$Wilcoxon_p < 0.1 & results_all$Wilcoxon_p >= 0.05, na.rm = TRUE)))
    cat(sprintf("  Средняя точность RF: %.3f\n", mean(results_all$RF_Accuracy, na.rm = TRUE)))
    
    # ---- РАСПРЕДЕЛЕНИЕ ПО ТИПАМ ----
    cat("\n📊 Распределение по типам комбинаций:\n")
    type_stats <- table(results_all$Type)
    for(t in names(type_stats)) {
      cat(sprintf("  %s: %d комбинаций\n", t, type_stats[t]))
    }
    
    # Сохраняем
    write.table(results_all, "drug_combinations_results.tsv", 
                sep = "\t", row.names = FALSE, quote = FALSE)
    cat("\n✅ Результаты сохранены в: drug_combinations_results.tsv\n")
    
  } else {
    cat("\n❌ Нет результатов для анализа\n")
  }
  
  cat("\n=========================================")
  cat("\nАНАЛИЗ ЗАВЕРШЕН!")
  cat("\n=========================================\n")
  
  
  # ============================================================================
  # ГРАФИК: Fecal metabolome VS. Drug combinations
  # ============================================================================
  
  library(ggplot2)
  
  cat("\n--- Создание графика: drug_combined_plot.pdf ---\n")
  
  # --- Подготовка данных ---
  plot_data_drug <- results_all
  
  # Убираем NA
  plot_data_drug <- plot_data_drug[!is.na(plot_data_drug$Wilcoxon_p) & 
                                     !is.na(plot_data_drug$RF_Accuracy), ]
  
  # Создаем цветовые группы
  plot_data_drug$color_group <- ifelse(plot_data_drug$Wilcoxon_p < 0.05, "#00BFC4", "gray70")
  plot_data_drug$color_group <- ifelse(plot_data_drug$Wilcoxon_p < 0.05 & 
                                         plot_data_drug$RF_Accuracy > 0.85, "#F8766D", 
                                       plot_data_drug$color_group)
  
  # Добавляем подписи
  plot_data_drug$label <- ifelse(plot_data_drug$Wilcoxon_p < 0.05 | 
                                   plot_data_drug$RF_Accuracy > 0.85, 
                                 plot_data_drug$Feature, "")
  
  # --- СОЗДАНИЕ PDF ---
  pdf("Desktop/WORK/gut/1_stage/R/drug_combined_plot.pdf", width = 8, height = 5.5)
  
  p <- ggplot(plot_data_drug, aes(x = Wilcoxon_p, 
                                  y = RF_Accuracy,
                                  size = N_Positive + N_Negative,
                                  color = color_group)) +
    geom_point(alpha = 0.8) +
    geom_vline(xintercept = 0.05, linetype = "dashed", color = "red", size = 0.8) +
    geom_hline(yintercept = 0.85, linetype = "dashed", color = "red", size = 0.8) +
    geom_text(aes(label = label), 
              vjust = -0.8, 
              hjust = 0.5, 
              size = 3.5,
              check_overlap = TRUE) +
    labs(title = "Fecal metabolome VS. Drug combinations",
         x = "P-value",
         y = "Accuracy (Random Forest)") +
    theme_minimal() +
    theme(
      axis.text = element_text(size = 17),
      plot.title = element_text(hjust = 0.5, size = 17, face = "bold"),
      axis.text.y = element_text(size = 14),
      axis.text.x = element_text(size = 14),
      legend.position = "none"
    ) +
    xlim(0, 1) +
    scale_color_identity()
  
  print(p)
  dev.off()
  
  cat("  ✅ drug_combined_plot.pdf\n")
  
  