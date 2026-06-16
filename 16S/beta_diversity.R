#Не используем, берем второй вариант с усреднением



library(vegan)
library(ggplot2)
library(dplyr)

# 1. Data loading
data_rel <- read.table("Desktop/WORK/gut/1_stage/R/input.tsv", header = TRUE, sep = "\t", check.names = FALSE, fill = TRUE)
metadata <- read.table("Desktop/WORK/gut/1_stage/R/metadata.tsv", header = TRUE, sep = "\t", check.names = FALSE, fill = TRUE)

# 2. Extract Sample IDs and find common intersection
ids_rel  <- as.character(data_rel[, 1])
ids_meta <- as.character(metadata[, 1])
common_samples <- intersect(ids_meta, ids_rel)

if (length(common_samples) == 0) {
  stop("Error: No matching Sample IDs found!")
}

# 3. Positional alignment
metadata_aligned <- metadata[match(common_samples, ids_meta), , drop = FALSE]
data_rel_aligned <- data_rel[match(common_samples, ids_rel), , drop = FALSE]

# Isolate numeric relative abundance matrix
mat_rel <- as.matrix(sapply(data_rel_aligned[, -1], as.numeric))

# 4. Build a completely fresh, sterile clean data frame for analysis
model_df <- data.frame(
  diagnosis = as.character(metadata_aligned$diagnosis),
  patient   = as.character(metadata_aligned$patient),
  stringsAsFactors = FALSE
)

# Convert to clean factors with defined levels
model_df$diagnosis <- factor(model_df$diagnosis, levels = c("CD", "UC", "Control"))
model_df$patient   <- as.factor(model_df$patient)

# Filter out any rows with NA values in the key metadata fields
valid_rows <- !is.na(model_df$diagnosis) & !is.na(model_df$patient)
model_df_cleaned <- model_df[valid_rows, ]
mat_rel_cleaned  <- mat_rel[valid_rows, ]

# 5. Calculate Bray-Curtis distance matrix
bray_dist <- vegdist(mat_rel_cleaned, method = "bray")

# 6. PERMANOVA with repeated measures control
set.seed(42)
permanova_res <- adonis2(bray_dist ~ diagnosis, data = model_df_cleaned, strata = model_df_cleaned$patient)
print(permanova_res)

# Extract PERMANOVA metrics safely
p_permanova  <- head(na.omit(permanova_res[["Pr(>F)"]]), 1)
r2_permanova <- head(na.omit(permanova_res[["R2"]]), 1)

# 7. Homogeneity of Multivariate Dispersions (PERMDISP)
permdisp_res <- betadisper(bray_dist, group = model_df_cleaned$diagnosis)

set.seed(42)
permdisp_test <- permutest(permdisp_res, permutations = 999)
print(permdisp_test)

# Extract GLOBAL PERMDISP p-value
p_permdisp_global <- head(na.omit(permdisp_test$tab[["Pr(>F)"]]), 1)

# --- НОВЫЙ БЛОК: ПОПАРНЫЙ ТЕСТ PERMDISP С КОРРЕКЦИЕЙ ХОЛМА ---
# Извлекаем расстояния от точек до центроидов своих групп
distances_to_centroid <- permdisp_res$distances

# Проводим попарные Т-тесты (так как это линейные расстояния)
pairwise_disp <- pairwise.t.test(distances_to_centroid, model_df_cleaned$diagnosis, 
                                 p.adjust.method = "holm")
print(pairwise_disp)

# Извлекаем попарные скорректированные p-values из матрицы результатов
p_matrix <- pairwise_disp$p.value
p_disp_CD_UC      <- p_matrix["UC", "CD"]
p_disp_CD_Control <- p_matrix["Control", "CD"]
p_disp_UC_Control <- p_matrix["Control", "UC"]

# Функция для красивого форматирования p-values
fmt_p <- function(p) {
  if (is.na(p)) return("NA")
  if (p < 0.001) return(sprintf("%.2e", p))
  return(as.character(round(p, 3)))
}
# -------------------------------------------------------------

# 8. Principal Coordinate Analysis (PCoA / CmdScale)
pcoa_res <- cmdscale(bray_dist, k = 2, eig = TRUE)

# Calculate variance percentages
all_eigenvalues <- pcoa_res$eig
positive_eigenvalues <- all_eigenvalues[all_eigenvalues > 0]
total_positive_sum <- sum(positive_eigenvalues)

# Calculate axes variance metrics safely
pc_percentages <- (positive_eigenvalues / total_positive_sum) * 100
pc1_lab <- sprintf("PCoA 1 (%.1f%%)", head(pc_percentages, 1))
pc2_lab <- sprintf("PCoA 2 (%.1f%%)", head(tail(pc_percentages, -1), 1))

# Build a plotting data frame
pcoa_plot_df <- data.frame(
  PC1       = pcoa_res$points[, 1],
  PC2       = pcoa_res$points[, 2],
  diagnosis = model_df_cleaned$diagnosis,
  patient   = model_df_cleaned$patient
)

# Форматируем новый расширенный блок текста для подписи графика
caption_text <- sprintf(
  "PERMANOVA: R² = %.3f, p-value = %.3f\nPERMDISP: p-value = %.3f\nPairwise p-value adj.: CD vs UC = %s, CD vs Control = %s, UC vs Control = %s", 
  r2_permanova, 
  p_permanova, 
  p_permdisp_global,
  fmt_p(p_disp_CD_UC),
  fmt_p(p_disp_CD_Control),
  fmt_p(p_disp_UC_Control)
)

# 9. Data Visualization (PCoA Ordination Plot with Advanced Bottom Stats)
p_pcoa <- ggplot(pcoa_plot_df, aes(x = PC1, y = PC2, color = diagnosis, fill = diagnosis)) +
  stat_ellipse(geom = "polygon", alpha = 0.1, aes(group = diagnosis), linewidth = 0.5) +
  geom_point(size = 3, alpha = 0.8) +
  scale_color_manual(values = c("CD" = "#D62828", "UC" = "#7209B7", "Control" = "#2A9D8F")) +
  scale_fill_manual(values = c("CD" = "#D62828", "UC" = "#7209B7", "Control" = "#2A9D8F")) +
  labs(
    x = pc1_lab,
    y = pc2_lab,
    title = "Beta Diversity (Bray-Curtis)",
    caption = caption_text,
    color = "Diagnosis",
    fill = "Diagnosis"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title   = element_text(size = 16, face = "bold", hjust = 0.5),
    plot.caption = element_text(size = 10.5, face = "bold.italic", hjust = 0.5, color = "grey20", margin = margin(t = 18)),
    axis.title   = element_text(size = 16),
    axis.text    = element_text(size = 14),
    legend.position = "right",
    panel.grid.major = element_line(linewidth = 0.2, colour = "grey90"),
    panel.grid.minor = element_line(linewidth = 0.12, colour = "grey93")
  )

# 10. Save PCoA plot as PDF
pdf("Desktop/WORK/gut/1_stage/R/beta_diversity.pdf", width = 8, height = 5.5) # Немного расширили холст под длинную подпись
print(p_pcoa)
dev.off()



















library(vegan)
library(ggplot2)
library(dplyr)

# 1. Загрузка данных
data_rel <- read.table("Desktop/WORK/gut/1_stage/R/input.tsv", header = TRUE, sep = "\t", check.names = FALSE, fill = TRUE)
metadata <- read.table("Desktop/WORK/gut/1_stage/R/metadata.tsv", header = TRUE, sep = "\t", check.names = FALSE, fill = TRUE)

# 2. Выравнивание образцов
ids_rel  <- as.character(data_rel[, 1])
ids_meta <- as.character(metadata[, 1])
common_samples <- intersect(ids_meta, ids_rel)

if (length(common_samples) == 0) {
  stop("Ошибка: совпадений Sample ID не найдено!")
}

metadata_aligned <- metadata[match(common_samples, ids_meta), , drop = FALSE]
data_rel_aligned <- data_rel[match(common_samples, ids_rel), , drop = FALSE]

# 3. Объединение метаданных и обилий бактерий
full_df <- data.frame(
  patient   = as.character(metadata_aligned$patient),
  diagnosis = as.character(metadata_aligned$diagnosis),
  stringsAsFactors = FALSE
)

# Переводим обилие бактерий в численный формат
mat_numeric <- sapply(data_rel_aligned[, -1], as.numeric)
full_df <- cbind(full_df, mat_numeric)

# 4. Усреднение (коллапсирование) профилей по каждому пациенту
collapsed_df <- full_df %>%
  filter(!is.na(diagnosis) & !is.na(patient)) %>%
  group_by(patient, diagnosis) %>%
  summarise(across(everything(), mean, na.rm = TRUE), .groups = 'drop')

# Превращаем диагноз в фактор с нужным порядком групп
collapsed_df$diagnosis <- factor(collapsed_df$diagnosis, levels = c("CD", "UC", "Control"))

# Разделяем обратно на чистые метаданные и матрицу обилий
model_df_cleaned <- collapsed_df[, c("patient", "diagnosis")]
mat_rel_cleaned  <- as.matrix(collapsed_df[, -c(1, 2)])

# 5. Расчет матрицы расстояний Брея-Кёртиса
bray_dist <- vegdist(mat_rel_cleaned, method = "bray")

# 6. Глобальный PERMANOVA (теперь без ошибок, так как повторов больше нет)
set.seed(42)
permanova_res <- adonis2(bray_dist ~ diagnosis, data = model_df_cleaned, permutations = 999)
print("--- Глобальный PERMANOVA ---")
print(permanova_res)

p_permanova  <- head(na.omit(permanova_res[["Pr(>F)"]]), 1)
r2_permanova <- head(na.omit(permanova_res[["R2"]]), 1)


# 7. ПОПАРНЫЙ ТЕСТ PERMANOVA (Пост-хок с коррекцией Холма)
pairwise_permanova <- function(distance_mat, meta_df, g1, g2) {
  # FIX: Changed '%in=' to '%in%'
  keep_idx <- meta_df$diagnosis %in% c(g1, g2) 
  sub_meta <- meta_df[keep_idx, ]
  sub_meta$diagnosis <- droplevels(factor(sub_meta$diagnosis))
  sub_dist <- as.dist(as.matrix(distance_mat)[keep_idx, keep_idx])
  
  res <- adonis2(sub_dist ~ diagnosis, data = sub_meta, permutations = 999)
  return(head(na.omit(res[["Pr(>F)"]]), 1))
}


set.seed(42)
p_raw_cd_uc      <- pairwise_permanova(bray_dist, model_df_cleaned, "CD", "UC")
p_raw_cd_control <- pairwise_permanova(bray_dist, model_df_cleaned, "CD", "Control")
p_raw_uc_control <- pairwise_permanova(bray_dist, model_df_cleaned, "UC", "Control")

# Применяем коррекцию Холма на множественные сравнения
p_perm_adj <- p.adjust(c(p_raw_cd_uc, p_raw_cd_control, p_raw_uc_control), method = "holm")
p_perm_CD_UC      <- p_perm_adj[1]
p_perm_CD_Control <- p_perm_adj[2]
p_perm_UC_Control <- p_perm_adj[3]

# 8. Анализ однородности дисперсий (PERMDISP)
permdisp_res <- betadisper(bray_dist, group = model_df_cleaned$diagnosis)

set.seed(42)
permdisp_test <- permutest(permdisp_res, permutations = 999)
print("--- Глобальный PERMDISP ---")
print(permdisp_test)

p_permdisp_global <- head(na.omit(permdisp_test$tab[["Pr(>F)"]]), 1)

# Попарный PERMDISP тест через расстояния до центроидов
distances_to_centroid <- permdisp_res$distances
pairwise_disp <- pairwise.t.test(distances_to_centroid, model_df_cleaned$diagnosis, p.adjust.method = "holm")
print("--- Попарный PERMDISP ---")
print(pairwise_disp)

p_disp_matrix <- pairwise_disp$p.value
p_disp_CD_UC      <- p_disp_matrix["UC", "CD"]
p_disp_CD_Control <- p_disp_matrix["Control", "CD"]
p_disp_UC_Control <- p_disp_matrix["Control", "UC"]

# Функция для форматирования p-values на графике
fmt_p <- function(p) {
  if (is.na(p)) return("NA")
  if (p < 0.001) return(sprintf("%.2e", p))
  return(as.character(round(p, 3)))
}

# Строка подписи для графика
caption_text <- sprintf(
  "PERMANOVA: R² = %.3f, p-value = %.3f\nPERMDISP: p-value = %.3f (CD vs UC=%s, CD vs Ctrl=%s, UC vs Ctrl=%s)\nPairwise PERMANOVA adj.: CD vs UC = %s, CD vs Control = %s, UC vs Control = %s", 
  r2_permanova, p_permanova, p_permdisp_global,
  fmt_p(p_disp_CD_UC), fmt_p(p_disp_CD_Control), fmt_p(p_disp_UC_Control),
  fmt_p(p_perm_CD_UC), fmt_p(p_perm_CD_Control), fmt_p(p_perm_UC_Control)
)

# 9. Построение оординации PCoA
pcoa_res <- cmdscale(bray_dist, k = 2, eig = TRUE)
total_positive_sum <- sum(pcoa_res$eig[pcoa_res$eig > 0])
pc_percentages <- (pcoa_res$eig[pcoa_res$eig > 0] / total_positive_sum) * 100

pc1_lab <- sprintf("PCoA 1 (%.1f%%)", pc_percentages[1])
pc2_lab <- sprintf("PCoA 2 (%.1f%%)", pc_percentages[2])

pcoa_plot_df <- data.frame(
  PC1       = pcoa_res$points[, 1],
  PC2       = pcoa_res$points[, 2],
  diagnosis = model_df_cleaned$diagnosis
)

# 10. Визуализация и сохранение в PDF
p_pcoa <- ggplot(pcoa_plot_df, aes(x = PC1, y = PC2, color = diagnosis, fill = diagnosis)) +
  stat_ellipse(geom = "polygon", alpha = 0.1, aes(group = diagnosis), linewidth = 0.5) +
  geom_point(size = 3, alpha = 0.8) +
  scale_color_manual(values = c("CD" = "#D62828", "UC" = "#7209B7", "Control" = "#2A9D8F")) +
  scale_fill_manual(values = c("CD" = "#D62828", "UC" = "#7209B7", "Control" = "#2A9D8F")) +
  labs(
    x = pc1_lab, y = pc2_lab,
    title = "Beta Diversity (Bray-Curtis)",
    caption = caption_text,
    color = "Diagnosis", fill = "Diagnosis"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title   = element_text(size = 16, face = "bold", hjust = 0.5),
    plot.caption = element_text(size = 9.5, face = "bold.italic", hjust = 0.5, color = "grey20", margin = margin(t = 18)),
    axis.title   = element_text(size = 16),
    axis.text    = element_text(size = 14),
    legend.position = "right",
    panel.grid.major = element_line(linewidth = 0.2, colour = "grey90")
  )

pdf("Desktop/WORK/gut/1_stage/R/beta_diversity.pdf", width = 8.5, height = 6)
print(p_pcoa)
dev.off()
