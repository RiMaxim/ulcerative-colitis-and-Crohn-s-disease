# Suppress warnings from lme4/reformulas
options(warn = -1)  # Temporarily disable warnings

library(maaslin3)
library(dplyr)
library(ggplot2)
library(forcats)
library(tidyr)
library(patchwork)

# ==============================================================================
# PART 1: RUNNING MAASLIN3
# ==============================================================================

# 1. Load data
input_data <- read.table("Desktop/github/input.tsv", header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)
input_metadata <- read.table("Desktop/github/metadata.tsv", header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)

# 2. Create binary variable (CD vs Control)
input_metadata$CD_vs_Control <- ifelse(
  input_metadata$diagnosis == "CD", "CD",
  ifelse(input_metadata$diagnosis == "Control", "Control", NA)
)

# Remove rows with NA (exclude UC patients)
metadata_filtered <- input_metadata[!is.na(input_metadata$CD_vs_Control), ]
metadata_filtered$patient <- as.character(metadata_filtered$patient)

# Align rows
common_samples <- intersect(rownames(metadata_filtered), rownames(input_data))
metadata_filtered <- metadata_filtered[common_samples, , drop = FALSE]
data_filtered <- input_data[common_samples, , drop = FALSE]

# 3. Run Maaslin3 with fixes
cat("Running Maaslin3 analysis...\n")
fit_data <- maaslin3(
  input_data     = data_filtered,
  input_metadata = metadata_filtered,
  output         = "Desktop/github/maaslin3_results_CD",
  fixed_effects  = c("CD_vs_Control", "age", "sex"),
  random_effects = c("patient"),
  reference      = c("CD_vs_Control=Control", "sex=F"), 
  normalization  = "TSS",
  transform      = "LOG",
  correction     = "BH",                            
  min_abundance  = 0.0001,
  min_prevalence = 0.1,
  small_random_effects = TRUE
)

# Re-enable warnings for the rest of the code
options(warn = 0)

# 4. Load generated results table
results <- read.delim("Desktop/github/maaslin3_results_CD/all_results.tsv", check.names = FALSE)

# Filter significant features by joint q-value
res_sig <- results %>% filter(pval_joint < 0.05)
res_cd  <- res_sig %>% filter(metadata == "CD_vs_Control")

p_abund <- NULL
p_prev  <- NULL

# LEFT PANEL: ABUNDANCE MODEL
abund <- res_cd %>%
  filter(model == "abundance") %>% 
  # Remove rows with critical NAs
  filter(!is.na(coef), !is.na(stderr), !is.na(pval_joint)) %>%
  mutate(
    feature = fct_reorder(feature, coef),
    direction = ifelse(coef > 0, "Higher in CD", "Lower in CD"),
    lower = coef - 1.96 * stderr,
    upper = coef + 1.96 * stderr,
    neglogp = -log10(pval_joint)
  )

if (nrow(abund) > 0) {
  p_abund <- ggplot(abund, aes(x = coef, y = feature, color = direction)) +
    geom_point(aes(size = neglogp), alpha = 0.8) +
    geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.2, alpha = 0.6) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    scale_color_manual(values = c("Higher in CD" = "#1B4F72", "Lower in CD" = "#C06014")) +
    scale_size_continuous(name = "-log10(p-value)") +
    labs(
      title = "Abundance associations",
      subtitle = "CD vs Control, adjusted for age, sex, and repeated measures",
      x = "β coefficient (±95% CI)", y = ""
    ) +
    theme_minimal(base_family = "Helvetica") +
    theme(
      plot.title.position = "plot",   
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),   
      plot.subtitle = element_text(size = 14, hjust = 0.5, color = "grey30"),              
      legend.position = "right",
      legend.direction = "vertical",   
      legend.justification = "top",    
      axis.text.y = element_text(face = "italic", size = 12),
      axis.title.x = element_text(size = 16),
      axis.text.x = element_text(size = 14),
      legend.text = element_text(size = 14),
      legend.title = element_text(size = 16)
    )
}

# RIGHT PANEL: PREVALENCE MODEL
prev <- res_cd %>%
  filter(model == "prevalence") %>% 
  mutate(
    feature = fct_reorder(feature, coef),
    direction = ifelse(coef > 0, "More frequent in CD", "Less frequent in CD"),
    lower = coef - 1.96 * stderr,
    upper = coef + 1.96 * stderr,
    neglogp = -log10(pval_joint)
  )

if (nrow(prev) > 0) {
  p_prev <- ggplot(prev, aes(x = coef, y = feature, color = direction)) +
    geom_point(aes(size = neglogp), alpha = 0.8) +
    geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.2, alpha = 0.6) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    scale_color_manual(values = c("More frequent in CD" = "#2E8B57", "Less frequent in CD" = "#B22222")) +
    scale_size_continuous(name = "-log10(p-value)") +
    labs(
      title = "Prevalence associations",
      subtitle = "CD vs Control, adjusted for age, sex, and repeated measures",
      x = "β coefficient (±95% CI)", y = ""
    ) +
    theme_minimal(base_family = "Helvetica") +
    theme(
      plot.title.position = "plot",   
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),   
      plot.subtitle = element_text(size = 14, hjust = 0.5, color = "grey30"),               
      legend.position = "right",
      legend.direction = "vertical",   
      legend.justification = "top",    
      axis.text.y = element_text(face = "italic", size = 12),
      axis.title.x = element_text(size = 16),
      axis.text.x = element_text(size = 14),
      legend.text = element_text(size = 14),
      legend.title = element_text(size = 16)
    )
}

# COMBINE AND SAVE THE MULTI-PANEL PLOT
if (!is.null(p_abund) && !is.null(p_prev)) {
  
  combined_dge <- (p_abund | p_prev) + 
    plot_annotation(
      title = "Bacterial biomarkers associated with Crohn's disease",
      theme = theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5, margin = margin(b = 15)))
    )
  
  ggsave("Desktop/github/CD_biomarkers.pdf", combined_dge, device = cairo_pdf, width =16, height = 10)
  cat("Success: Multi-panel plot saved to Desktop/github/CD_biomarkers.pdf\n")
  
} else if (!is.null(p_abund) || !is.null(p_prev)) {
  
  only_plot <- if(is.null(p_abund)) p_prev else p_abund
  ggsave("Desktop/github/CD_biomarkers.pdf", only_plot, device = cairo_pdf, width = 8, height = 10)
  cat("Only one model had significant results. Single panel saved to Desktop/github/CD_biomarkers.pdf\n")
  
} else {
  cat("No significant features found under p-value < 0.05. PDF was not generated.\n")
}












# Suppress warnings from lme4/reformulas
options(warn = -1)  # Temporarily disable warnings

library(maaslin3)
library(dplyr)
library(ggplot2)
library(forcats)
library(tidyr)
library(patchwork)

# ==============================================================================
# PART 1: RUNNING MAASLIN3
# ==============================================================================

# 1. Load data
input_data <- read.table("Desktop/github/input.tsv", header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)
input_metadata <- read.table("Desktop/github/metadata.tsv", header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)

# 2. Create binary variable (UC vs Control)
input_metadata$UC_vs_Control <- ifelse(
  input_metadata$diagnosis == "UC", "UC",
  ifelse(input_metadata$diagnosis == "Control", "Control", NA)
)

# Remove rows with NA (exclude CD patients)
metadata_filtered <- input_metadata[!is.na(input_metadata$UC_vs_Control), ]
metadata_filtered$patient <- as.character(metadata_filtered$patient)

# Align rows
common_samples <- intersect(rownames(metadata_filtered), rownames(input_data))
metadata_filtered <- metadata_filtered[common_samples, , drop = FALSE]
data_filtered <- input_data[common_samples, , drop = FALSE]

# 3. Run Maaslin3 with fixes
cat("Running Maaslin3 analysis...\n")
fit_data <- maaslin3(
  input_data     = data_filtered,
  input_metadata = metadata_filtered,
  output         = "Desktop/github/maaslin3_results_UC",
  fixed_effects  = c("UC_vs_Control", "age", "sex"),
  random_effects = c("patient"),
  reference      = c("UC_vs_Control=Control", "sex=F"), 
  normalization  = "TSS",
  transform      = "LOG",
  correction     = "BH",                            
  min_abundance  = 0.0001,
  min_prevalence = 0.1,
  small_random_effects = TRUE
)

# Re-enable warnings for the rest of the code
options(warn = 0)

# 4. Load generated results table
results <- read.delim("Desktop/github/maaslin3_results_UC/all_results.tsv", check.names = FALSE)

# Filter significant features by joint q-value
res_sig <- results %>% filter(pval_joint < 0.05)
res_uc  <- res_sig %>% filter(metadata == "UC_vs_Control")

p_abund <- NULL
p_prev  <- NULL

# LEFT PANEL: ABUNDANCE MODEL
abund <- res_uc %>%
  filter(model == "abundance") %>% 
  # Remove rows with critical NAs
  filter(!is.na(coef), !is.na(stderr), !is.na(pval_joint)) %>%
  mutate(
    feature = fct_reorder(feature, coef),
    direction = ifelse(coef > 0, "Higher in UC", "Lower in UC"),
    lower = coef - 1.96 * stderr,
    upper = coef + 1.96 * stderr,
    neglogp = -log10(pval_joint)
  )

if (nrow(abund) > 0) {
  p_abund <- ggplot(abund, aes(x = coef, y = feature, color = direction)) +
    geom_point(aes(size = neglogp), alpha = 0.8) +
    geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.2, alpha = 0.6) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    scale_color_manual(values = c("Higher in UC" = "#1B4F72", "Lower in UC" = "#C06014")) +
    scale_size_continuous(name = "-log10(p-value)") +
    labs(
      title = "Abundance associations",
      subtitle = "UC vs Control, adjusted for age, sex, and repeated measures",
      x = "β coefficient (±95% CI)", y = ""
    ) +
    theme_minimal(base_family = "Helvetica") +
    theme(
      plot.title.position = "plot",   
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),   
      plot.subtitle = element_text(size = 14, hjust = 0.5, color = "grey30"),              
      legend.position = "right",
      legend.direction = "vertical",   
      legend.justification = "top",    
      axis.text.y = element_text(face = "italic", size = 12),
      axis.title.x = element_text(size = 16),
      axis.text.x = element_text(size = 14),
      legend.text = element_text(size = 14),
      legend.title = element_text(size = 16)
    )
}

# RIGHT PANEL: PREVALENCE MODEL
prev <- res_uc %>%
  filter(model == "prevalence") %>% 
  mutate(
    feature = fct_reorder(feature, coef),
    direction = ifelse(coef > 0, "More frequent in UC", "Less frequent in UC"),
    lower = coef - 1.96 * stderr,
    upper = coef + 1.96 * stderr,
    neglogp = -log10(pval_joint)
  )

if (nrow(prev) > 0) {
  p_prev <- ggplot(prev, aes(x = coef, y = feature, color = direction)) +
    geom_point(aes(size = neglogp), alpha = 0.8) +
    geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.2, alpha = 0.6) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    scale_color_manual(values = c("More frequent in UC" = "#2E8B57", "Less frequent in UC" = "#B22222")) +
    scale_size_continuous(name = "-log10(p-value)") +
    labs(
      title = "Prevalence associations",
      subtitle = "UC vs Control, adjusted for age, sex, and repeated measures",
      x = "β coefficient (±95% CI)", y = ""
    ) +
    theme_minimal(base_family = "Helvetica") +
    theme(
      plot.title.position = "plot",   
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),   
      plot.subtitle = element_text(size = 14, hjust = 0.5, color = "grey30"),              
      legend.position = "right",
      legend.direction = "vertical",   
      legend.justification = "top",    
      axis.text.y = element_text(face = "italic", size = 12),
      axis.title.x = element_text(size = 16),
      axis.text.x = element_text(size = 14),
      legend.text = element_text(size = 14),
      legend.title = element_text(size = 16)
    )
}

# COMBINE AND SAVE THE MULTI-PANEL PLOT
if (!is.null(p_abund) && !is.null(p_prev)) {
  
  combined_dge <- (p_abund | p_prev) + 
    plot_annotation(
      title = "Bacterial biomarkers associated with Ulcerative Colitis",
      theme = theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5, margin = margin(b = 15)))
    )
  
  ggsave("Desktop/github/UC_biomarkers.pdf", combined_dge, device = cairo_pdf, width = 16, height = 10)
  cat("Success: Multi-panel plot saved to Desktop/github/UC_biomarkers.pdf\n")
  
} else if (!is.null(p_abund) || !is.null(p_prev)) {
  
  only_plot <- if(is.null(p_abund)) p_prev else p_abund
  ggsave("Desktop/github/UC_biomarkers.pdf", only_plot, device = cairo_pdf, width = 8, height = 10)
  cat("Only one model had significant results. Single panel saved to Desktop/github/UC_biomarkers.pdf\n")
  
} else {
  cat("No significant features found under p-value < 0.05. PDF was not generated.\n")
}


