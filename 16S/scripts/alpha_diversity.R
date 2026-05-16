library(lme4)
library(lmerTest)
library(emmeans)
library(vegan)
library(ggplot2)
library(ggsignif)
library(dplyr)
library(patchwork)

# 1. Data loading
data_rel    <- read.table("Desktop/github/input.tsv", header = TRUE, sep = "\t", check.names = FALSE, fill = TRUE)
data_counts <- read.table("Desktop/github/input_counts.tsv", header = TRUE, sep = "\t", check.names = FALSE, fill = TRUE)
metadata    <- read.table("Desktop/github/metadata.tsv", header = TRUE, sep = "\t", check.names = FALSE, fill = TRUE)

# 2. Extract Sample IDs and find common intersection
ids_rel    <- as.character(data_rel[, 1])
ids_counts <- as.character(data_counts[, 1])
ids_meta   <- as.character(metadata[, 1])

common_samples <- intersect(intersect(ids_meta, ids_rel), ids_counts)

if (length(common_samples) == 0) {
  stop("Error: No matching Sample IDs found across all input files!")
}

# 3. Positional alignment
metadata_aligned <- metadata[match(common_samples, ids_meta), , drop = FALSE]
data_rel_aligned <- data_rel[match(common_samples, ids_rel), , drop = FALSE]
data_cnt_aligned <- data_counts[match(common_samples, ids_counts), , drop = FALSE]

mat_rel <- as.matrix(sapply(data_rel_aligned[, -1], as.numeric))
mat_cnt <- as.matrix(sapply(data_cnt_aligned[, -1], as.numeric))
mat_cnt_int <- round(mat_cnt, 0)

# 4. Alpha diversity calculations
shannon_vec  <- diversity(mat_rel, index = "shannon")
richness_vec <- specnumber(mat_cnt_int)

# 5. Build clean data frame
model_df <- data.frame(
  shannon   = as.numeric(shannon_vec),
  richness  = as.numeric(richness_vec),
  diagnosis = as.character(metadata_aligned$diagnosis),
  patient   = as.character(metadata_aligned$patient),
  stringsAsFactors = FALSE
)
model_df$diagnosis <- factor(model_df$diagnosis, levels = c("CD", "UC", "Control"))
model_df <- na.omit(model_df)

# --- AUTOMATED MODELING AND PLOTTING FUNCTION ---
run_analysis_and_plot <- function(metric_name, y_label, fill_colors) {
  
  # Fit Linear Mixed-Effects Model
  formula_str <- paste(metric_name, "~ diagnosis + (1 | patient)")
  mix_model   <- lmer(as.formula(formula_str), data = model_df)
  
  # Pairwise comparisons WITH Holm correction
  emm          <- emmeans(mix_model, ~ diagnosis)
  lmm_summary  <- as.data.frame(pairs(emm, adjust = "holm"))
  
  # Pairwise comparisons WITHOUT correction (Raw P-values)
  lmm_raw      <- as.data.frame(pairs(emm, adjust = "none"))
  
  # Helper to extract values (fixed to fetch BOTH raw and adjusted correctly)
  get_p_vals <- function(g1, g2) {
    row_adj <- lmm_summary[grepl(g1, lmm_summary$contrast) & grepl(g2, lmm_summary$contrast), ]
    row_raw <- lmm_raw[grepl(g1, lmm_raw$contrast) & grepl(g2, lmm_raw$contrast), ]
    return(list(raw = row_raw$p.value, adj = row_adj$p.value))
  }
  
  p_CD_UC      <- get_p_vals("CD", "UC")
  p_UC_Control <- get_p_vals("UC", "Control")
  p_CD_Control <- get_p_vals("CD", "Control")
  
  # Format double label function (Raw P over Adjusted P)
  fmt_double_p <- function(p_list) {
    p_raw <- p_list$raw
    p_adj <- p_list$adj
    
    if (is.na(p_raw) || is.na(p_adj)) return("NA")
    
    txt_raw <- if (p_raw < 0.001) sprintf("%.2e", p_raw) else round(p_raw, 2)
    txt_adj <- if (p_adj < 0.001) sprintf("%.2e", p_adj) else round(p_adj, 2)
    
    return(paste0("p-value = ", txt_raw, "\np-value adj. = ", txt_adj))
  }
  
  # Ordered lists from bottom bracket to top bracket
  p_labels    <- c(fmt_double_p(p_CD_UC), fmt_double_p(p_UC_Control), fmt_double_p(p_CD_Control))
  comparisons <- list(c("CD", "UC"), c("UC", "Control"), c("CD", "Control"))
  
  y_max       <- max(model_df[[metric_name]], na.rm = TRUE)
  y_min       <- min(model_df[[metric_name]], na.rm = TRUE)
  
  # Step size for brackets
  step        <- (y_max - y_min) * 0.18 
  y_positions <- y_max + c(step, step * 2.3, step * 3.6)
  
  # Generate Plot
  plt <- ggplot(model_df, aes(x = diagnosis, y = .data[[metric_name]], fill = diagnosis)) +
    geom_violin(trim = FALSE, alpha = 0.6, width = 0.8) +
    geom_boxplot(width = 0.15, fill = "white", outlier.shape = NA, alpha = 0.8) +
    scale_fill_manual(values = fill_colors) +
    labs(title = "Alpha Diversity", x = "", y = y_label) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title   = element_text(size = 16, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 16),
      axis.text  = element_text(size = 14),
      legend.position = "none",
      panel.grid.major = element_line(linewidth = 0.2, colour = "grey90"),
      panel.grid.minor = element_line(linewidth = 0.12, colour = "grey93")
    ) +
    geom_signif(comparisons = comparisons,
                annotations = p_labels,
                y_position = y_positions,
                tip_length = 0.02,
                textsize = 3.0,       
                vjust = -0.1) +
    # Adds dynamic breathing room above the top bracket so labels do not get cut off
    expand_limits(y = max(y_positions) + step)
  
  return(plt)
}

# 6. Generate both plots
colors <- c("CD" = "#D62828", "UC" = "#7209B7", "Control" = "#2A9D8F")
plot_shannon  <- run_analysis_and_plot("shannon", "Shannon Index", colors)
plot_richness <- run_analysis_and_plot("richness", "Observed Species (Richness)", colors)

# 7. Combine plots side-by-side and save
combined_plot <- plot_shannon + plot_richness

pdf("Desktop/github/alpha_diversity.pdf", width = 8, height = 5.5) 
print(combined_plot)
dev.off()

