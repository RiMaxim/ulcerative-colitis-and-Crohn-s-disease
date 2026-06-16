library(maaslin3)
library(dplyr)
library(ggplot2)
library(forcats)
library(tidyr)
library(patchwork)

# ==============================================================================
# PART 1: LOAD DATA
# ==============================================================================

input_data <- read.table(
  "Desktop/WORK/gut/1_stage/R/input.tsv",
  header = TRUE, row.names = 1,
  sep = "\t", check.names = FALSE
)

input_metadata <- read.table(
  "Desktop/WORK/gut/1_stage/R/metadata.tsv",
  header = TRUE, row.names = 1,
  sep = "\t", check.names = FALSE
)

# ==============================================================================
# PART 2: CREATE BINARY GROUP (CD vs Control)
# ==============================================================================

input_metadata$CD_vs_Control <- ifelse(
  input_metadata$diagnosis == "CD", "CD",
  ifelse(input_metadata$diagnosis == "Control", "Control", NA)
)

# remove UC
metadata_filtered <- input_metadata[!is.na(input_metadata$CD_vs_Control), ]

# basic cleaning
metadata_filtered$patient <- as.character(metadata_filtered$patient)
metadata_filtered$age <- as.numeric(as.character(metadata_filtered$age))

# ==============================================================================
# PART 3: RENAME DRUG VARIABLES (IMPORTANT FIX)
# ==============================================================================
colnames(metadata_filtered)[colnames(metadata_filtered) == "drug_1"] <- "asa5"
colnames(metadata_filtered)[colnames(metadata_filtered) == "drug_2"] <- "gc"
colnames(metadata_filtered)[colnames(metadata_filtered) == "drug_3"] <- "thiopurines"
colnames(metadata_filtered)[colnames(metadata_filtered) == "drug_4"] <- "biologics"
colnames(metadata_filtered)[colnames(metadata_filtered) == "drug_5"] <- "jak"
colnames(metadata_filtered)[colnames(metadata_filtered) == "drug_6"] <- "antibiotics"

# ensure numeric (important for MaAsLin3)
drug_cols <- c("asa5","gc","thiopurines","biologics","jak","antibiotics")
metadata_filtered[drug_cols] <- lapply(metadata_filtered[drug_cols], function(x) as.numeric(as.character(x)))

lapply(metadata_filtered[drug_cols], table, useNA = "ifany")

# ==============================================================================
# PART 4: ALIGN DATA
# ==============================================================================

common_samples <- intersect(rownames(metadata_filtered), rownames(input_data))

metadata_filtered <- metadata_filtered[common_samples, , drop = FALSE]
data_filtered <- input_data[common_samples, , drop = FALSE]

metadata_filtered$CD_vs_Control <- factor(metadata_filtered$CD_vs_Control, levels = c("Control", "CD"))
metadata_filtered$sex <- factor(metadata_filtered$sex, levels = c("F", "M"))

# ==============================================================================
# PART 5: RUN MAASLIN3 (UPDATED MODEL)
# ==============================================================================

cat("Running MaAsLin3 analysis...\n")

fit_data <- maaslin3(
  input_data     = data_filtered,
  input_metadata = metadata_filtered,
  output         = "Desktop/WORK/gut/1_stage/R/maaslin3_results_CD",
  
  fixed_effects  = c(
    "CD_vs_Control",
    "age",
    "sex",
    
    # DRUGS (UPDATED)
    "asa5",
    "gc",
    "thiopurines",
    "biologics",
    "jak",
    "antibiotics"
  ),
  
  random_effects = c("patient"),
  
  reference = c(
    "CD_vs_Control=Control",
    "sex=F"
  ),
  
  normalization = "TSS",
  transform     = "LOG",
  correction    = "BH",
  
  min_abundance  = 0.0001,
  min_prevalence = 0.1,
  
  small_random_effects = TRUE
)

# ==============================================================================
# PART 6: LOAD RESULTS
# ==============================================================================

results <- read.delim(
  "Desktop/WORK/gut/1_stage/R/maaslin3_results_CD/all_results.tsv",
  check.names = FALSE
)

# keep significant
res_sig <- results %>% filter(pval_joint < 0.05)
res_cd  <- res_sig %>% filter(metadata == "CD_vs_Control")

# ==============================================================================
# PART 7: ABUNDANCE PLOT
# ==============================================================================

abund <- res_cd %>%
  filter(model == "abundance") %>%
  filter(!is.na(coef), !is.na(stderr), !is.na(pval_joint)) %>%
  mutate(
    feature = fct_reorder(feature, coef),
    direction = ifelse(coef > 0, "Higher in CD", "Lower in CD"),
    lower = coef - 1.96 * stderr,
    upper = coef + 1.96 * stderr,
    neglogp = -log10(pval_joint)
  )

p_abund <- ggplot(abund, aes(x = coef, y = feature, color = direction)) +
  geom_point(aes(size = neglogp), alpha = 0.8) +
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.2, alpha = 0.6) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_color_manual(values = c("Higher in CD" = "#1B4F72", "Lower in CD" = "#C06014")) +
  scale_size_continuous(name = "-log10(p-value)") +
  labs(
    title = "Abundance associations (CD vs Control)",
    subtitle = "CD vs Control, adjusted for age, sex, repeated measures and 6 drug classes",
    x = "β coefficient (±95% CI)", y = ""
  ) +
  theme_minimal()

# ==============================================================================
# PART 8: PREVALENCE PLOT
# ==============================================================================

prev <- res_cd %>%
  filter(model == "prevalence") %>%
  mutate(
    feature = fct_reorder(feature, coef),
    direction = ifelse(coef > 0, "More frequent in CD", "Less frequent in CD"),
    lower = coef - 1.96 * stderr,
    upper = coef + 1.96 * stderr,
    neglogp = -log10(pval_joint)
  )

p_prev <- ggplot(prev, aes(x = coef, y = feature, color = direction)) +
  geom_point(aes(size = neglogp), alpha = 0.8) +
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.2, alpha = 0.6) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_color_manual(values = c("More frequent in CD" = "#2E8B57", "Less frequent in CD" = "#B22222")) +
  scale_size_continuous(name = "-log10(p-value)") +
  labs(
    title = "Prevalence associations",
    subtitle = "CD vs Control, adjusted for age, sex, repeated measures and 6 drug classes",
    x = "β coefficient (±95% CI)", y = ""
  ) +
  theme_minimal()

# ==============================================================================
# PART 9: COMBINE & SAVE
# ==============================================================================

combined <- (p_abund | p_prev) +
  plot_annotation(
    title = "Bacterial biomarkers associated with Crohn's disease",
    theme = theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5, margin = margin(b = 15)))
  )

ggsave(
  "Desktop/WORK/gut/1_stage/R/CD_biomarkers_new.pdf",
  combined,
  width = 16, height = 10, device = cairo_pdf
)

cat("DONE: MaAsLin3 analysis with drug adjustment completed\n")










library(maaslin3)
library(dplyr)
library(ggplot2)
library(forcats)
library(tidyr)
library(patchwork)

# ==============================================================================
# PART 1: LOAD DATA
# ==============================================================================

input_data <- read.table(
  "Desktop/WORK/gut/1_stage/R/input.tsv",
  header = TRUE, row.names = 1,
  sep = "\t", check.names = FALSE
)

input_metadata <- read.table(
  "Desktop/WORK/gut/1_stage/R/metadata.tsv",
  header = TRUE, row.names = 1,
  sep = "\t", check.names = FALSE
)

# ==============================================================================
# PART 2: CREATE UC vs CONTROL
# ==============================================================================

input_metadata$UC_vs_Control <- ifelse(
  input_metadata$diagnosis == "UC", "UC",
  ifelse(input_metadata$diagnosis == "Control", "Control", NA)
)

# keep only UC + Control
metadata_filtered <- input_metadata[!is.na(input_metadata$UC_vs_Control), ]

# basic cleaning
metadata_filtered$patient <- as.character(metadata_filtered$patient)
metadata_filtered$age <- as.numeric(as.character(metadata_filtered$age))

# ==============================================================================
# PART 3: RENAME DRUG VARIABLES
# ==============================================================================

colnames(metadata_filtered)[colnames(metadata_filtered) == "drug_1"] <- "asa5"
colnames(metadata_filtered)[colnames(metadata_filtered) == "drug_2"] <- "gc"
colnames(metadata_filtered)[colnames(metadata_filtered) == "drug_3"] <- "thiopurines"
colnames(metadata_filtered)[colnames(metadata_filtered) == "drug_4"] <- "biologics"
colnames(metadata_filtered)[colnames(metadata_filtered) == "drug_5"] <- "jak"
colnames(metadata_filtered)[colnames(metadata_filtered) == "drug_6"] <- "antibiotics"

drug_cols <- c("asa5","gc","thiopurines","biologics","jak","antibiotics")

metadata_filtered[drug_cols] <- lapply(metadata_filtered[drug_cols], function(x) {
  x <- as.numeric(x)
  x[is.na(x)] <- 0
  x
})

# ==============================================================================
# PART 4: ALIGN DATA
# ==============================================================================

common_samples <- intersect(rownames(metadata_filtered), rownames(input_data))

metadata_filtered <- metadata_filtered[common_samples, , drop = FALSE]
data_filtered <- input_data[common_samples, , drop = FALSE]

metadata_filtered$UC_vs_Control <- factor(metadata_filtered$UC_vs_Control, levels = c("Control", "UC"))
metadata_filtered$sex <- factor(metadata_filtered$sex, levels = c("F", "M"))

# ==============================================================================
# PART 5: RUN MAASLIN3
# ==============================================================================

cat("Running MaAsLin3 UC analysis...\n")

fit_data <- maaslin3(
  input_data     = data_filtered,
  input_metadata = metadata_filtered,
  output         = "Desktop/WORK/gut/1_stage/R/maaslin3_results_UC",
  
  fixed_effects  = c(
    "UC_vs_Control",
    "age",
    "sex",
    
    # drugs
    "asa5",
    "gc",
    "thiopurines",
    "biologics",
    "jak",
    "antibiotics"
  ),
  
  random_effects = c("patient"),
  
  reference = c(
    "UC_vs_Control=Control",
    "sex=F"
  ),
  
  normalization = "TSS",
  transform     = "LOG",
  correction    = "BH",
  
  min_abundance  = 0.0001,
  min_prevalence = 0.1,
  
  small_random_effects = TRUE
)

# ==============================================================================
# PART 6: LOAD RESULTS
# ==============================================================================

results <- read.delim(
  "Desktop/WORK/gut/1_stage/R/maaslin3_results_UC/all_results.tsv",
  check.names = FALSE
)

res_sig <- results %>% filter(pval_joint < 0.05)
res_uc  <- res_sig %>% filter(metadata == "UC_vs_Control")

# ==============================================================================
# PART 7: PLOT (ABUNDANCE)
# ==============================================================================

abund <- res_uc %>%
  filter(model == "abundance") %>%
  filter(!is.na(coef), !is.na(stderr), !is.na(pval_joint)) %>%
  mutate(
    feature = fct_reorder(feature, coef),
    direction = ifelse(coef > 0, "Higher in UC", "Lower in UC"),
    lower = coef - 1.96 * stderr,
    upper = coef + 1.96 * stderr,
    neglogp = -log10(pval_joint)
  )

p_abund <- ggplot(abund, aes(x = coef, y = feature, color = direction)) +
  geom_point(aes(size = neglogp), alpha = 0.8) +
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.2, alpha = 0.6) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_color_manual(values = c("Higher in UC" = "#1B4F72", "Lower in UC" = "#C06014")) +
  scale_size_continuous(name = "-log10(p-value)") +
  labs(
    title = "Abundance associations",
    subtitle = "UC vs Control, adjusted for age, sex, repeated measures and 6 drug classes",
    x = "β coefficient (±95% CI)", y = ""
  ) +
  theme_minimal()

# ==============================================================================
# PART 8: PLOT (PREVALENCE)
# ==============================================================================

prev <- res_uc %>%
  filter(model == "prevalence") %>%
  mutate(
    feature = fct_reorder(feature, coef),
    direction = ifelse(coef > 0, "More frequent in UC", "Less frequent in UC"),
    lower = coef - 1.96 * stderr,
    upper = coef + 1.96 * stderr,
    neglogp = -log10(pval_joint)
  )

p_prev <- ggplot(prev, aes(x = coef, y = feature, color = direction)) +
  geom_point(aes(size = neglogp), alpha = 0.8) +
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.2, alpha = 0.6) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_color_manual(values = c("More frequent in UC" = "#2E8B57", "Less frequent in UC" = "#B22222")) +
  scale_size_continuous(name = "-log10(p-value)") +
  labs(
    title = "Prevalence associations",
    subtitle = "UC vs Control, adjusted for age, sex, repeated measures and 6 drug classes",
    x = "β coefficient (±95% CI)", y = ""
  ) +
  theme_minimal()

# ==============================================================================
# PART 9: SAVE
# ==============================================================================

combined <- (p_abund | p_prev) +
  plot_annotation(
    title = "Bacterial biomarkers associated with Ulcerative Colitis",
    theme = theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5, margin = margin(b = 15)))
  )


ggsave(
  "Desktop/WORK/gut/1_stage/R/UC_biomarkers_new.pdf",
  combined,
  width = 16, height = 10, device = cairo_pdf
)

cat("DONE: UC MaAsLin3 analysis completed\n")




