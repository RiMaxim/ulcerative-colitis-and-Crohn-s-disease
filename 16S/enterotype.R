library(vegan)
library(cluster)
library(ggplot2)
library(dplyr)
library(patchwork)
# 1. Fixed random generator for total reproducibility
set.seed(42)

# 2. Data loading
data_rel <- read.table("Desktop/WORK/gut/1_stage/R/input.tsv", header = TRUE, sep = "\t", check.names = FALSE, fill = TRUE)
metadata <- read.table("Desktop/WORK/gut/1_stage/R/metadata.tsv", header = TRUE, sep = "\t", check.names = FALSE, fill = TRUE)

# 3. Sample alignment
ids_rel  <- as.character(data_rel[, 1])
ids_meta <- as.character(metadata[, 1])
common_samples <- intersect(ids_meta, ids_rel)

if (length(common_samples) == 0) {
  stop("Error: No matching Sample IDs found!")
}

metadata_aligned <- metadata[match(common_samples, ids_meta), , drop = FALSE]
data_rel_aligned <- data_rel[match(common_samples, ids_rel), , drop = FALSE]

# Isolate clean relative abundance matrix
mat_rel <- as.matrix(sapply(data_rel_aligned[, -1], as.numeric))

# Filter out low-quality rows if any exist in factors
valid_rows <- !is.na(metadata_aligned$diagnosis) & !is.na(metadata_aligned$patient)
metadata_cleaned <- metadata_aligned[valid_rows, ]
mat_rel_cleaned  <- mat_rel[valid_rows, ]

# 4. Custom function to calculate Jensen-Shannon Divergence (JSD) Matrix
# JSD is the golden standard metric for gut microbiota enterotyping
dist_jsd <- function(matrix_data) {
  # Small epsilon to avoid logarithm of zero
  matrix_data <- matrix_data + 1e-10
  matrix_data <- matrix_data / rowSums(matrix_data)
  
  n_samples <- nrow(matrix_data)
  jsd_matrix <- matrix(0, nrow = n_samples, ncol = n_samples)
  rownames(jsd_matrix) <- rownames(matrix_data)
  colnames(jsd_matrix) <- rownames(matrix_data)
  
  for (i in 1:(n_samples - 1)) {
    for (j in (i + 1):n_samples) {
      p <- matrix_data[i, ]
      q <- matrix_data[j, ]
      m <- 0.5 * (p + q)
      
      kl_pm <- sum(p * log2(p / m))
      kl_qm <- sum(q * log2(q / m))
      
      jsd_val <- sqrt(0.5 * kl_pm + 0.5 * kl_qm)
      jsd_matrix[i, j] <- jsd_val
      jsd_matrix[j, i] <- jsd_val
    }
  }
  return(as.dist(jsd_matrix))
}

jsd_distances <- dist_jsd(mat_rel_cleaned)

# 5. Determine optimal number of clusters (K) using Silhouette Index
# We test partition layouts from 2 to 6 clusters
max_k <- 6
silhouette_scores <- numeric(max_k)

for (k in 2:max_k) {
  pam_fit <- pam(jsd_distances, diss = TRUE, k = k)
  silhouette_scores[k] <- pam_fit$silinfo$avg.width
}

# Auto-select K that maximizes the Silhouette Width
optimal_k <- which.max(silhouette_scores)
cat(sprintf("Optimal number of Enterotypes detected: %d (Silhouette: %.3f)\n", 
            optimal_k, silhouette_scores[optimal_k]))

# 6. Final clustering assignment using the optimal K layout
final_pam <- pam(jsd_distances, diss = TRUE, k = optimal_k)
metadata_cleaned$enterotype <- paste0("ET_", final_pam$clustering)

# 7. Coordinate extraction for plotting (PCoA based on JSD)
pcoa_jsd <- cmdscale(jsd_distances, k = 2, eig = TRUE)
pcoa_eig <- pcoa_jsd$eig
pcoa_eig_pos <- pcoa_eig[pcoa_eig > 0]
total_sum <- sum(pcoa_eig_pos)

pc1_label <- sprintf("PCoA 1 (%.1f%%)", (head(pcoa_eig_pos, 1) / total_sum) * 100)
pc2_label <- sprintf("PCoA 2 (%.1f%%)", (head(tail(pcoa_eig_pos, -1), 1) / total_sum) * 100)

plot_df <- data.frame(
  PC1        = pcoa_jsd$points[, 1],
  PC2        = pcoa_jsd$points[, 2],
  enterotype = metadata_cleaned$enterotype,
  diagnosis  = factor(metadata_cleaned$diagnosis, levels = c("CD", "UC", "Control"))
)

# 8. Plot 1: Ordination space colored by detected Enterotypes
p_et <- ggplot(plot_df, aes(x = PC1, y = PC2, color = enterotype, fill = enterotype)) +
  stat_ellipse(geom = "polygon", alpha = 0.1, linewidth = 0.5) +
  geom_point(size = 3, alpha = 0.8) +
  scale_color_brewer(palette = "Set1") +
  scale_fill_brewer(palette = "Set1") +
  labs(x = pc1_label, y = pc2_label, title = "Detected Gut Enterotypes (JSD Matrix)", color = "Enterotype", fill = "Enterotype") +
  theme_minimal(base_size = 14) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        axis.title   = element_text(size = 16),
        axis.text    = element_text(size = 14))

# 9. Plot 2: Distribution of Enterotypes across clinical diagnoses
p_dist <- ggplot(plot_df, aes(x = diagnosis, fill = enterotype)) +
  geom_bar(position = "fill", width = 0.6) +
  scale_fill_brewer(palette = "Set1") +
  labs(x = "", y = "Proportion of Samples", title = "Enterotype Distribution by Group", fill = "Enterotype") +
  theme_minimal(base_size = 14) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
  axis.title   = element_text(size = 16),
  axis.text    = element_text(size = 14))
# 10. Save plots as PDF
pdf("Desktop/WORK/gut/1_stage/R/enterotypes_1.pdf", width = 8, height = 5.5)
print(p_et)
dev.off()
pdf("Desktop/WORK/gut/1_stage/R/enterotypes_2.pdf", width = 8, height = 5.5)
p_dist
dev.off()
# 11. Perform Fisher's exact test for 3x3 contingency table
fisher_res <- fisher.test(table(metadata_cleaned$diagnosis, metadata_cleaned$enterotype), simulate.p.value = TRUE, B = 10000)
print(fisher_res)

# Print p-value to console
cat(sprintf("\nStatistical significance of enterotype distribution (Fisher's Test): P = %.4f\n", fisher_res$p.value))








library(dplyr)
library(tidyr)
library(ggplot2)
library(ggalluvial) # Если не установлен, запустите: install.packages("ggalluvial")

# 1. Сборка лонгитюдного датафрейма (берем только пациентов с повторными визитами)
trajectory_df <- metadata_cleaned %>%
  select(patient, diagnosis, enterotype) %>%
  filter(!is.na(diagnosis) & !is.na(patient) & enterotype %in% c("ET_1", "ET_2")) %>%
  group_by(patient) %>%
  # Фильтруем: оставляем только тех, кто пришел минимум 2 раза
  filter(n() >= 2) %>%
  # Автоматически нумеруем визиты хронологически по порядку строк
  mutate(visit_order = paste0("Visit_", row_number())) %>%
  ungroup()

# Переводим визиты в фактор, чтобы зафиксировать их порядок на графике
trajectory_df$visit_order <- factor(trajectory_df$visit_order, 
                                    levels = paste0("Visit_", 1:max(as.numeric(gsub("Visit_", "", trajectory_df$visit_order)))))
trajectory_df$diagnosis <- factor(trajectory_df$diagnosis, levels = c("Control", "CD", "UC"))

# 2. Расчет парных переходов (Текущий шаг -> Следующий шаг)
transition_pairs <- trajectory_df %>%
  group_by(patient, diagnosis) %>%
  mutate(
    current_et = enterotype,
    next_et    = lead(enterotype)
  ) %>%
  filter(!is.na(next_et)) %>% # Исключаем последний визит, у которого нет пары
  ungroup()

cat(sprintf("Всего успешно зафиксировано переходов во времени: %d\n", nrow(transition_pairs)))

# 3. Функция генерации Матрицы Переходов (%) для конкретной группы
get_transition_matrix <- function(df, group_name) {
  sub_df <- df %>% filter(diagnosis == group_name)
  
  if(nrow(sub_df) == 0) {
    return(matrix(NA, nrow=2, ncol=2, dimnames=list(c("ET_1","ET_2"), c("ET_1","ET_2"))))
  }
  
  # Строим таблицу переходов и нормируем по строкам (margin = 1) для получения вероятностей
  raw_table <- table(sub_df$current_et, sub_df$next_et)
  prob_table <- prop.table(raw_table, margin = 1) * 100
  return(round(prob_table, 1))
}

# Выводим матрицы в консоль
cat("\n=== МАТРИЦА ПЕРЕХОДОВ ДЛЯ ГРУППЫ CONTROL (%) ===\n")
print(get_transition_matrix(transition_pairs, "Control"))

cat("\n=== МАТРИЦА ПЕРЕХОДОВ ДЛЯ ГРУППЫ CD (%) ===\n")
print(get_transition_matrix(transition_pairs, "CD"))

cat("\n=== МАТРИЦА ПЕРЕХОДОВ ДЛЯ ГРУППЫ UC (%) ===\n")
print(get_transition_matrix(transition_pairs, "UC"))


# 4. Визуализация: Отрисовка траекторий (Alluvial Diagram) по группам
# Ограничим график первыми тремя визитами, чтобы он оставался читаемым
plot_alluvial_df <- trajectory_df %>%
  filter(visit_order %in% c("Visit_1", "Visit_2", "Visit_3")) %>%
  group_by(diagnosis, visit_order, enterotype) %>%
  summarise(Freq = n(), .groups = 'drop')

p_alluvial <- ggplot(plot_alluvial_df,
                     aes(x = visit_order, stratum = enterotype, alluvium = enterotype, y = Freq, fill = enterotype)) +
  geom_flow(alpha = 0.4, color = "darkgrey", linewidth = 0.2) +
  geom_stratum(alpha = 0.9, width = 0.4, color = "black") +
  facet_wrap(~diagnosis, scales = "free_y") +
  scale_fill_manual(values = c("ET_1" = "#D62828", "ET_2" = "blue")) +
  labs(
    title = "Longitudinal Dynamics of Enterotypes",
    subtitle = "Tracking microbial community shifts across sequential patient visits",
    x = "Timeline",
    y = "Number of Samples",
    fill = "Enterotype"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title    = element_text(face = "bold", hjust = 0.5, size = 16),
    plot.subtitle = element_text(hjust = 0.5, size = 12, face = "italic"),
    strip.text    = element_text(face = "bold", size = 14, color = "grey10"),
    legend.position = "bottom",
    panel.grid.major.x = element_blank()
  ) +ylim(0,30)

# Сохраняем аллювиальный график в PDF
pdf("Desktop/WORK/gut/1_stage/R/enterotype_3.pdf", width = 8, height = 5.5)
print(p_alluvial)
dev.off()






