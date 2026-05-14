library(vegan)
library(cluster)
library(ggplot2)
library(dplyr)
library(patchwork)
# 1. Fixed random generator for total reproducibility
set.seed(42)

# 2. Data loading
data_rel <- read.table("Desktop/github/input.tsv", header = TRUE, sep = "\t", check.names = FALSE, fill = TRUE)
metadata <- read.table("Desktop/github/metadata.tsv", header = TRUE, sep = "\t", check.names = FALSE, fill = TRUE)

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
pdf("Desktop/github/enterotypes_1.pdf", width = 8, height = 5.5)
print(p_et)
dev.off()
pdf("Desktop/github/enterotypes_2.pdf", width = 8, height = 5.5)
p_dist
dev.off()
# 11. Perform Fisher's exact test for 3x3 contingency table
fisher_res <- fisher.test(table(metadata_cleaned$diagnosis, metadata_cleaned$enterotype), simulate.p.value = TRUE, B = 10000)
print(fisher_res)

# Print p-value to console
cat(sprintf("\nStatistical significance of enterotype distribution (Fisher's Test): P = %.4f\n", fisher_res$p.value))
