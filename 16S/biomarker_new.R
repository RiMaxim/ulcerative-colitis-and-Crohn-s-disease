library(maaslin3)
library(dplyr)
library(ggplot2)
library(forcats)
library(tidyr)
library(patchwork)

input_data <- read.table("Desktop/WORK/gut/1_stage/R/input.tsv", header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)
input_metadata <- read.table("Desktop/WORK/gut/1_stage/R/metadata.tsv", header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)

dim(input_data)
dim(input_metadata)







###########################################################################PART1

# ==============================================================================
# PART 2: CREATE BINARY GROUP (CD vs Control)
# ==============================================================================

input_metadata$CD_vs_Control <- ifelse(
  input_metadata$diagnosis == "CD", "CD",
  ifelse(input_metadata$diagnosis == "Control", "Control", NA)
)

# remove UC
metadata_filtered <- input_metadata[!is.na(input_metadata$CD_vs_Control), ]
dim(metadata_filtered)
# ==============================================================================
# PART 4: ALIGN DATA
# ==============================================================================

common_samples <- intersect(rownames(metadata_filtered), rownames(input_data))

metadata_filtered <- metadata_filtered[common_samples, , drop = FALSE]
data_filtered <- input_data[common_samples, , drop = FALSE]
dim(data_filtered)

metadata_filtered$CD_vs_Control <- factor(metadata_filtered$CD_vs_Control, levels = c("Control", "CD"))
metadata_filtered$sex <- factor(metadata_filtered$sex, levels = c("F", "M"))
metadata_filtered$patient <- as.character(metadata_filtered$patient)
metadata_filtered$age <- as.numeric(as.character(metadata_filtered$age))

# ==============================================================================
# PART 5: RUN MAASLIN3
# ==============================================================================

cat("Running MaAsLin3 analysis...\n")

fit_data <- maaslin3(
  input_data     = data_filtered,
  input_metadata = metadata_filtered,
  output         = "Desktop/WORK/gut/1_stage/R/maaslin3_results_CD",
  
  fixed_effects  = c(
    "CD_vs_Control",
    "age",
    "sex"),
  
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
    direction = ifelse(coef > 0, "Higher in CD", "Higher in Control"),
    lower = coef - 1.96 * stderr,
    upper = coef + 1.96 * stderr,
    neglogp = -log10(pval_joint)
  )

p_abund <- ggplot(abund, aes(x = coef, y = feature, color = direction)) +
  geom_point(aes(size = neglogp), alpha = 0.8) +
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.2, alpha = 0.6) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_color_manual(
  name = "Direction",
  values = c("Higher in CD" = "#1B4F72", "Higher in Control" = "#C06014")
) +
  scale_size_continuous(
    name = "-log10(P-value)",
  ) +
  scale_y_discrete(
    labels = function(x) gsub("_", " ", x)
  ) +
  labs(
    title = "Abundance associations",
    subtitle = "CD vs Control\nadjusted for age and gender",
    x = "β coefficient (±95% CI)", y = ""
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    legend.box = "vertical",
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    axis.text.y = element_text(face = "italic")
  )+
  guides(
    color = guide_legend(order = 1),
    size = guide_legend(order = 2)
  )

# ==============================================================================
# PART 8: PREVALENCE PLOT
# ==============================================================================

prev <- res_cd %>%
  filter(model == "prevalence") %>%
  mutate(
    feature = fct_reorder(feature, coef),
    direction = ifelse(coef > 0, "More frequent in CD", "More frequent in Control"),
    lower = coef - 1.96 * stderr,
    upper = coef + 1.96 * stderr,
    neglogp = -log10(pval_joint)
  )

p_prev <- ggplot(prev, aes(x = coef, y = feature, color = direction)) +
  geom_point(aes(size = neglogp), alpha = 0.8) +
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.2, alpha = 0.6) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_color_manual(
  name = "Direction",
  values = c("More frequent in CD" = "#2E8B57", "More frequent in Control" = "#B22222")
) +
  scale_size_continuous(
    name = "-log10(P-value)",
  ) +
  scale_y_discrete(
    labels = function(x) gsub("_", " ", x)
  ) +
  labs(
    title = "Prevalence associations",
    subtitle = "CD vs Control\nadjusted for age and gender",
    x = "β coefficient (±95% CI)", y = ""
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    legend.box = "vertical",
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    axis.text.y = element_text(face = "italic")
  )+
  guides(
    color = guide_legend(order = 1),
    size = guide_legend(order = 2)
  )

# ==============================================================================
# PART 9: COMBINE & SAVE
# ==============================================================================

combined1 <- (p_abund | p_prev) +
  plot_annotation(
    title = "Bacterial biomarkers associated with Crohn's disease",
    theme = theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5, margin = margin(b = 15)))
  )















##########################################################################PART2

# ==============================================================================
# PART 2: CREATE UC vs CONTROL
# ==============================================================================

input_metadata$UC_vs_Control <- ifelse(
  input_metadata$diagnosis == "UC", "UC",
  ifelse(input_metadata$diagnosis == "Control", "Control", NA)
)

# keep only UC + Control
metadata_filtered <- input_metadata[!is.na(input_metadata$UC_vs_Control), ]

# ==============================================================================
# PART 4: ALIGN DATA
# ==============================================================================

common_samples <- intersect(rownames(metadata_filtered), rownames(input_data))

metadata_filtered <- metadata_filtered[common_samples, , drop = FALSE]
data_filtered <- input_data[common_samples, , drop = FALSE]

metadata_filtered$UC_vs_Control <- factor(metadata_filtered$UC_vs_Control, levels = c("Control", "UC"))
metadata_filtered$sex <- factor(metadata_filtered$sex, levels = c("F", "M"))
metadata_filtered$patient <- as.character(metadata_filtered$patient)
metadata_filtered$age <- as.numeric(as.character(metadata_filtered$age))


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
    "sex"),
  
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
    direction = ifelse(coef > 0, "Higher in UC", "Higher in Control"),
    lower = coef - 1.96 * stderr,
    upper = coef + 1.96 * stderr,
    neglogp = -log10(pval_joint)
  )

p_abund <- ggplot(abund, aes(x = coef, y = feature, color = direction)) +
  geom_point(aes(size = neglogp), alpha = 0.8) +
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.2, alpha = 0.6) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_color_manual(
    name = "Direction",
    values = c("Higher in UC" = "#1B4F72", "Higher in Control" = "#C06014")
  ) +
  scale_size_continuous(
    name = "-log10(P-value)",
  ) +
  scale_y_discrete(
    labels = function(x) gsub("_", " ", x)
  ) +
  labs(
    title = "Abundance associations",
    subtitle = "UC vs Control\nadjusted for age and gender",
    x = "β coefficient (±95% CI)", y = ""
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    legend.box = "vertical",
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    axis.text.y = element_text(face = "italic")
  )+
  guides(
    color = guide_legend(order = 1),
    size = guide_legend(order = 2)
  )

# ==============================================================================
# PART 8: PLOT (PREVALENCE)
# ==============================================================================

prev <- res_uc %>%
  filter(model == "prevalence") %>%
  mutate(
    feature = fct_reorder(feature, coef),
    direction = ifelse(coef > 0, "More frequent in UC", "More frequent in Control"),
    lower = coef - 1.96 * stderr,
    upper = coef + 1.96 * stderr,
    neglogp = -log10(pval_joint)
  )

p_prev <- ggplot(prev, aes(x = coef, y = feature, color = direction)) +
  geom_point(aes(size = neglogp), alpha = 0.8) +
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.2, alpha = 0.6) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_color_manual(
    name = "Direction",
    values = c("More frequent in UC" = "#2E8B57", "More frequent in Control" = "#B22222")
  ) +
  scale_size_continuous(
    name = "-log10(P-value)",
  ) +
  scale_y_discrete(
    labels = function(x) gsub("_", " ", x)
  ) +
  labs(
    title = "Prevalence associations",
    subtitle = "UC vs Control\nadjusted for age and gender",
    x = "β coefficient (±95% CI)", y = ""
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    legend.box = "vertical",
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    axis.text.y = element_text(face = "italic")
  )+
  guides(
    color = guide_legend(order = 1),
    size = guide_legend(order = 2)
  )
# ==============================================================================
# PART 9: SAVE
# ==============================================================================

combined2 <- (p_abund | p_prev) +
  plot_annotation(
    title = "Bacterial biomarkers associated with Ulcerative Colitis",
    theme = theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5, margin = margin(b = 15)))
  )












###########################################################################PART3



# ==============================================================================
# PART 2: CREATE CD vs UC GROUPS (ONLY CD AND UC, NO CONTROLS)
# ==============================================================================

# Оставляем только CD и UC (исключаем Control)
metadata_filtered <- input_metadata[input_metadata$diagnosis %in% c("CD", "UC"), ]

# Создаем переменную CD_vs_UC
metadata_filtered$CD_vs_UC <- factor(
  metadata_filtered$diagnosis,
  levels = c("UC", "CD")  # UC как референс
)
dim(metadata_filtered)

# ==============================================================================
# PART 3: RENAME DRUG VARIABLES
# ==============================================================================

colnames(metadata_filtered)[colnames(metadata_filtered) == "drug_1"] <- "5-ASA"
colnames(metadata_filtered)[colnames(metadata_filtered) == "drug_2"] <- "Glucocorticoids"
colnames(metadata_filtered)[colnames(metadata_filtered) == "drug_3"] <- "Immunosuppressors"
colnames(metadata_filtered)[colnames(metadata_filtered) == "drug_4"] <- "Genetically_Engineered_Therapies"
colnames(metadata_filtered)[colnames(metadata_filtered) == "drug_5"] <- "Inhibitor_JAK_kinase"
colnames(metadata_filtered)[colnames(metadata_filtered) == "drug_6"] <- "Antibiotics"

# ensure numeric
drug_cols <- c("5-ASA", "Glucocorticoids", "Immunosuppressors", 
               "Genetically_Engineered_Therapies", "Inhibitor_JAK_kinase", "Antibiotics")
metadata_filtered[drug_cols] <- lapply(metadata_filtered[drug_cols], function(x) as.numeric(as.character(x)))

metadata_filtered <- metadata_filtered[complete.cases(metadata_filtered[drug_cols]), ]
dim(metadata_filtered)
# ==============================================================================
# PART 4: ALIGN DATA
# ==============================================================================

common_samples <- intersect(rownames(metadata_filtered), rownames(input_data))

metadata_filtered <- metadata_filtered[common_samples, , drop = FALSE]
data_filtered <- input_data[common_samples, , drop = FALSE]

metadata_filtered$sex <- factor(metadata_filtered$sex, levels = c("F", "M"))
metadata_filtered$em <- factor(metadata_filtered$em, levels = c(0, 1))
metadata_filtered$patient <- as.character(metadata_filtered$patient)
metadata_filtered$age <- as.numeric(as.character(metadata_filtered$age))

# ==============================================================================
# PART 5: RUN MAASLIN3 (CD vs UC with drugs)
# ==============================================================================

cat("Running MaAsLin3 CD vs UC analysis with drugs...\n")

fit_data <- maaslin3(
  input_data     = data_filtered,
  input_metadata = metadata_filtered,
  output         = "Desktop/WORK/gut/1_stage/R/maaslin3_results_CD_vs_UC",
  
  fixed_effects  = c(
    "CD_vs_UC",
    "age",
    "sex",
    "em",
    "5-ASA",
    "Glucocorticoids",
    "Immunosuppressors",
    "Genetically_Engineered_Therapies",
    "Inhibitor_JAK_kinase",
    "Antibiotics"
  ),
  
  random_effects = c("patient"),
  
  reference = c(
    "CD_vs_UC=UC",  # UC как референс
    "sex=F",
    "em=0"
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
  "Desktop/WORK/gut/1_stage/R/maaslin3_results_CD_vs_UC/all_results.tsv",
  check.names = FALSE
)

# keep significant
res_sig <- results %>% filter(pval_joint < 0.05)
res_cd_vs_uc <- res_sig %>% filter(metadata == "CD_vs_UC")

# ==============================================================================
# PART 7: PLOT (ABUNDANCE)
# ==============================================================================

abund <- res_cd_vs_uc %>%
  filter(model == "abundance") %>%
  filter(!is.na(coef), !is.na(stderr), !is.na(pval_joint)) %>%
  mutate(
    feature = fct_reorder(feature, coef),
    direction = ifelse(coef > 0, "Higher in CD", "Higher in UC"),
    lower = coef - 1.96 * stderr,
    upper = coef + 1.96 * stderr,
    neglogp = -log10(pval_joint)
  )

p_abund <- ggplot(abund, aes(x = coef, y = feature, color = direction)) +
  geom_point(aes(size = neglogp), alpha = 0.8) +
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.2, alpha = 0.6) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_color_manual(
    name = "Direction",
    values = c("Higher in CD" = "#1B4F72", "Higher in UC" = "#C06014")
  ) +
  scale_size_continuous(
    name = "-log10(P-value)"
  ) +
  scale_y_discrete(
    labels = function(x) gsub("_", " ", x)
  ) +
  labs(
    title = "Abundance associations",
    subtitle = "CD vs UC\nadjusted for age, gender, therapy and extraintestinal manifistation",
    x = "β coefficient (±95% CI)", 
    y = ""
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    legend.box = "vertical",
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    axis.text.y = element_text(face = "italic")
  )+
  guides(
    color = guide_legend(order = 1),
    size = guide_legend(order = 2)
  )

# ==============================================================================
# PART 8: PLOT (PREVALENCE)
# ==============================================================================

prev <- res_cd_vs_uc %>%
  filter(model == "prevalence") %>%
  mutate(
    feature = fct_reorder(feature, coef),
    direction = ifelse(coef > 0, "More frequent in CD", "More frequent in UC"),
    lower = coef - 1.96 * stderr,
    upper = coef + 1.96 * stderr,
    neglogp = -log10(pval_joint)
  )

p_prev <- ggplot(prev, aes(x = coef, y = feature, color = direction)) +
  geom_point(aes(size = neglogp), alpha = 0.8) +
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.2, alpha = 0.6) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_color_manual(
    name = "Direction",
    values = c("More frequent in CD" = "#2E8B57", "More frequent in UC" = "#B22222")
  ) +
  scale_size_continuous(
    name = "-log10(P-value)"
  ) +
  scale_y_discrete(
    labels = function(x) gsub("_", " ", x)
  ) +
  labs(
    title = "Prevalence associations",
    subtitle = "CD vs UC\nadjusted for age, gender, therapy and extraintestinal manifistation",
    x = "β coefficient (±95% CI)", 
    y = ""
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    legend.box = "vertical",
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    axis.text.y = element_text(face = "italic")
  )+
  guides(
    color = guide_legend(order = 1),
    size = guide_legend(order = 2)
  )

# ==============================================================================
# PART 9: COMBINE & SAVE
# ==============================================================================

combined3 <- (p_abund | p_prev) +
  plot_annotation(
    title = "Bacterial biomarkers differentiating Crohn's disease and Ulcerative Colitis",
    theme = theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5, margin = margin(b = 15))
    )
  )




#####################################PART4


# ==============================================================================
# PART 2: FILTER CD PATIENTS AND CREATE ER VARIABLE
# ==============================================================================

# Оставляем только CD
metadata_filtered <- input_metadata[input_metadata$diagnosis == "CD", ]

# ==============================================================================
# PART 3: RENAME DRUG VARIABLES
# ==============================================================================

colnames(metadata_filtered)[colnames(metadata_filtered) == "drug_1"] <- "5-ASA"
colnames(metadata_filtered)[colnames(metadata_filtered) == "drug_2"] <- "Glucocorticoids"
colnames(metadata_filtered)[colnames(metadata_filtered) == "drug_3"] <- "Immunosuppressors"
colnames(metadata_filtered)[colnames(metadata_filtered) == "drug_4"] <- "Genetically_Engineered_Therapies"
colnames(metadata_filtered)[colnames(metadata_filtered) == "drug_5"] <- "Inhibitor_JAK_kinase"
colnames(metadata_filtered)[colnames(metadata_filtered) == "drug_6"] <- "Antibiotics"

# ensure numeric
drug_cols <- c("5-ASA", "Glucocorticoids", "Immunosuppressors", 
               "Genetically_Engineered_Therapies", "Inhibitor_JAK_kinase", "Antibiotics")
metadata_filtered[drug_cols] <- lapply(metadata_filtered[drug_cols], function(x) as.numeric(as.character(x)))
metadata_filtered <- metadata_filtered[complete.cases(metadata_filtered[drug_cols]), ]
dim(metadata_filtered)
# ==============================================================================
# PART 4: ALIGN DATA
# ==============================================================================

common_samples <- intersect(rownames(metadata_filtered), rownames(input_data))

metadata_filtered <- metadata_filtered[common_samples, , drop = FALSE]
data_filtered <- input_data[common_samples, , drop = FALSE]

metadata_filtered$sex <- factor(metadata_filtered$sex, levels = c("F", "M"))
metadata_filtered$er <- factor(metadata_filtered$er, levels = c(0, 1))
metadata_filtered$em <- factor(metadata_filtered$em, levels = c(0, 1))
metadata_filtered$patient <- as.character(metadata_filtered$patient)
metadata_filtered$age <- as.numeric(as.character(metadata_filtered$age))

# ==============================================================================
# PART 5: RUN MAASLIN3 (ER 1 vs 0 within CD, adjusting for EM)
# ==============================================================================

cat("Running MaAsLin3 CD-only ER analysis...\n")

fit_data <- maaslin3(
  input_data     = data_filtered,
  input_metadata = metadata_filtered,
  output         = "Desktop/WORK/gut/1_stage/R/maaslin3_results_CD_ER",
  
  fixed_effects  = c(
    "er",        # основная переменная интереса
    "em",        # вычитание/поправка на em
    "age",
    "sex",
    "5-ASA",
    "Glucocorticoids",
    "Immunosuppressors",
    "Genetically_Engineered_Therapies",
    "Inhibitor_JAK_kinase",
    "Antibiotics"
  ),
  
  random_effects = c("patient"),
  
  reference = c(
    "er=0",      # er=0 как референс
    "em=0",      # em=0 как референс
    "sex=F"
  ),
  
  normalization = "TSS",
  transform     = "LOG",
  correction    = "BH",
  
  min_abundance  = 0.0001,
  min_prevalence = 0.1,
  
  small_random_effects = TRUE
)

  #==============================================================================
  # PART 6: LOAD RESULTS
  # ==============================================================================

results <- read.delim(
  "Desktop/WORK/gut/1_stage/R/maaslin3_results_CD_ER/all_results.tsv",
  check.names = FALSE
)

# keep significant for er (main effect)
res_sig <- results %>% filter(pval_joint < 0.05)
res_cd_er <- res_sig %>% filter(metadata == "er")

# ==============================================================================
# PART 7: PLOT (ABUNDANCE) - ER 1 vs 0
# ==============================================================================

abund <- res_cd_er %>%
  filter(model == "abundance") %>%
  filter(!is.na(coef), !is.na(stderr), !is.na(pval_joint)) %>%
  mutate(
    feature = fct_reorder(feature, coef),
    direction = ifelse(coef > 0, "Higher in exacerbation", "Higher in remission"),
    lower = coef - 1.96 * stderr,
    upper = coef + 1.96 * stderr,
    neglogp = -log10(pval_joint)
  )

p_abund <- ggplot(abund, aes(x = coef, y = feature, color = direction)) +
  geom_point(aes(size = neglogp), alpha = 0.8) +
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.2, alpha = 0.6) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_color_manual(
    name = "Direction",
    values = c("Higher in exacerbation" = "#1B4F72", "Higher in remission" = "#C06014")
  ) +
  scale_size_continuous(
    name = "-log10(P-value)"
  ) +
  scale_y_discrete(
    labels = function(x) gsub("_", " ", x)
  ) +
  labs(
    title = "Abundance associations",
    subtitle = "CD patients only: exacerbation vs remission\nadjusted for age, gender, therapy and extraintestinal manifistation",
    x = "β coefficient (±95% CI)", 
    y = ""
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    legend.box = "vertical",
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    axis.text.y = element_text(face = "italic")
  )+
  guides(
    color = guide_legend(order = 1),
    size = guide_legend(order = 2)
  )

# ==============================================================================
# PART 8: PLOT (PREVALENCE)
# ==============================================================================

prev <- res_cd_er %>%
  filter(model == "prevalence") %>%
  mutate(
    feature = fct_reorder(feature, coef),
    direction = ifelse(coef > 0, "More frequent in exacerbation", "More frequent in remission"),
    lower = coef - 1.96 * stderr,
    upper = coef + 1.96 * stderr,
    neglogp = -log10(pval_joint)
  )

p_prev <- ggplot(prev, aes(x = coef, y = feature, color = direction)) +
  geom_point(aes(size = neglogp), alpha = 0.8) +
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.2, alpha = 0.6) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_color_manual(
    name = "Direction",
    values = c("More frequent in exacerbation" = "#2E8B57", "More frequent in remission" = "#B22222")
  ) +
  scale_size_continuous(
    name = "-log10(P-value)"
  ) +
  scale_y_discrete(
    labels = function(x) gsub("_", " ", x)
  ) +
  labs(
    title = "Prevalence associations",
    subtitle = "CD patients only: exacerbation vs remission\nadjusted for age, gender, therapy and extraintestinal manifistation",
    x = "β coefficient (±95% CI)", 
    y = ""
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    legend.box = "vertical",
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    axis.text.y = element_text(face = "italic")
  )+
  guides(
    color = guide_legend(order = 1),  # Direction ПЕРВАЯ
    size = guide_legend(order = 2)    # -log10(P-value) ВТОРАЯ
  )

# ==============================================================================
# PART 9: COMBINE & SAVE
# ==============================================================================

combined4 <- (p_abund | p_prev) +
  plot_annotation(
    title = "Bacterial biomarkers associated with exacerbation / remission status in Crohn's disease patients",
    theme = theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5, margin = margin(b = 15))
    )
  )



######################PART5


# ==============================================================================
# PART 2: FILTER UC PATIENTS AND CREATE ER VARIABLE
# ==============================================================================

# Оставляем только UC
metadata_filtered <- input_metadata[input_metadata$diagnosis == "UC", ]

# ==============================================================================
# PART 3: RENAME DRUG VARIABLES
# ==============================================================================

colnames(metadata_filtered)[colnames(metadata_filtered) == "drug_1"] <- "5-ASA"
colnames(metadata_filtered)[colnames(metadata_filtered) == "drug_2"] <- "Glucocorticoids"
colnames(metadata_filtered)[colnames(metadata_filtered) == "drug_3"] <- "Immunosuppressors"
colnames(metadata_filtered)[colnames(metadata_filtered) == "drug_4"] <- "Genetically_Engineered_Therapies"
colnames(metadata_filtered)[colnames(metadata_filtered) == "drug_5"] <- "Inhibitor_JAK_kinase"
colnames(metadata_filtered)[colnames(metadata_filtered) == "drug_6"] <- "Antibiotics"

# ensure numeric
drug_cols <- c("5-ASA", "Glucocorticoids", "Immunosuppressors", 
               "Genetically_Engineered_Therapies", "Inhibitor_JAK_kinase", "Antibiotics")
metadata_filtered[drug_cols] <- lapply(metadata_filtered[drug_cols], function(x) as.numeric(as.character(x)))
metadata_filtered <- metadata_filtered[complete.cases(metadata_filtered[drug_cols]), ]

# ==============================================================================
# PART 4: ALIGN DATA
# ==============================================================================

common_samples <- intersect(rownames(metadata_filtered), rownames(input_data))

metadata_filtered <- metadata_filtered[common_samples, , drop = FALSE]
data_filtered <- input_data[common_samples, , drop = FALSE]

metadata_filtered$sex <- factor(metadata_filtered$sex, levels = c("F", "M"))
metadata_filtered$er <- factor(metadata_filtered$er, levels = c(0, 1))
metadata_filtered$em <- factor(metadata_filtered$em, levels = c(0, 1))
metadata_filtered$patient <- as.character(metadata_filtered$patient)
metadata_filtered$age <- as.numeric(as.character(metadata_filtered$age))


# ==============================================================================
# PART 5: RUN MAASLIN3 (ER 1 vs 0 within UC, adjusting for EM)
# ==============================================================================

cat("Running MaAsLin3 UC-only ER analysis...\n")

fit_data <- maaslin3(
  input_data     = data_filtered,
  input_metadata = metadata_filtered,
  output         = "Desktop/WORK/gut/1_stage/R/maaslin3_results_UC_ER",
  
  fixed_effects  = c(
    "er",        # основная переменная интереса
    "em",        # вычитание/поправка на em
    "age",
    "sex",
    "5-ASA",
    "Glucocorticoids",
    "Immunosuppressors",
    "Genetically_Engineered_Therapies",
    "Inhibitor_JAK_kinase"
  ),
  
  random_effects = c("patient"),
  
  reference = c(
    "er=0",      # er=0 как референс
    "em=0",      # em=0 как референс
    "sex=F"
  ),
  
  normalization = "TSS",
  transform     = "LOG",
  correction    = "BH",
  
  min_abundance  = 0.0001,
  min_prevalence = 0.1,
  
  small_random_effects = TRUE
)



#==============================================================================
# PART 6: LOAD RESULTS
# ==============================================================================

results <- read.delim(
  "Desktop/WORK/gut/1_stage/R/maaslin3_results_UC_ER/all_results.tsv",
  check.names = FALSE
)

# keep significant for er (main effect)
res_sig <- results %>% filter(pval_joint < 0.05)
res_uc_er <- res_sig %>% filter(metadata == "er")

# ==============================================================================
# PART 7: PLOT (ABUNDANCE) - ER 1 vs 0
# ==============================================================================

abund <- res_uc_er %>%
  filter(model == "abundance") %>%
  filter(!is.na(coef), !is.na(stderr), !is.na(pval_joint)) %>%
  mutate(
    feature = fct_reorder(feature, coef),
    direction = ifelse(coef > 0, "Higher in exacerbation", "Higher in remission"),
    lower = coef - 1.96 * stderr,
    upper = coef + 1.96 * stderr,
    neglogp = -log10(pval_joint)
  )

p_abund <- ggplot(abund, aes(x = coef, y = feature, color = direction)) +
  geom_point(aes(size = neglogp), alpha = 0.8) +
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.2, alpha = 0.6) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_color_manual(
    name = "Direction",
    values = c("Higher in exacerbation" = "#1B4F72", "Higher in remission" = "#C06014")
  ) +
  scale_size_continuous(
    name = "-log10(P-value)"
  ) +
  scale_y_discrete(
    labels = function(x) gsub("_", " ", x)
  ) +
  labs(
    title = "Abundance associations",
    subtitle = "UC patients only: exacerbation vs remission\nadjusted for age, gender, therapy and extraintestinal manifistation",
    x = "β coefficient (±95% CI)", 
    y = ""
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    legend.box = "vertical",
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    axis.text.y = element_text(face = "italic")
  )+
  guides(
    color = guide_legend(order = 1),
    size = guide_legend(order = 2)
  )

# ==============================================================================
# PART 8: PLOT (PREVALENCE)
# ==============================================================================

prev <- res_uc_er %>%
  filter(model == "prevalence") %>%
  mutate(
    feature = fct_reorder(feature, coef),
    direction = ifelse(coef > 0, "More frequent in exacerbation", "More frequent in remission"),
    lower = coef - 1.96 * stderr,
    upper = coef + 1.96 * stderr,
    neglogp = -log10(pval_joint)
  )

p_prev <- ggplot(prev, aes(x = coef, y = feature, color = direction)) +
  geom_point(aes(size = neglogp), alpha = 0.8) +
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.2, alpha = 0.6) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_color_manual(
    name = "Direction",
    values = c("More frequent in exacerbation" = "#2E8B57", "More frequent in remission" = "#B22222")
  ) +
  scale_size_continuous(
    name = "-log10(P-value)"
  ) +
  scale_y_discrete(
    labels = function(x) gsub("_", " ", x)
  ) +
  labs(
    title = "Prevalence associations",
    subtitle = "UC patients only: exacerbation vs remission\nadjusted for age, gender, therapy and extraintestinal manifistation",
    x = "β coefficient (±95% CI)", 
    y = ""
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    legend.box = "vertical",
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    axis.text.y = element_text(face = "italic")
  )+
  guides(
    color = guide_legend(order = 1),
    size = guide_legend(order = 2)
  )

# ==============================================================================
# PART 9: COMBINE & SAVE
# ==============================================================================

combined5 <- (p_abund | p_prev) +
  plot_annotation(
    title = "Bacterial biomarkers associated with exacerbation / remission status in Ulcerative Colitis patients",
    theme = theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5, margin = margin(b = 15))
    )
  )


###################PART6




# ==============================================================================
# PART 2: FILTER CD PATIENTS AND CREATE ACTIVITY VARIABLE
# ==============================================================================

# Оставляем только CD
metadata_filtered <- input_metadata[input_metadata$diagnosis == "CD", ]

# Преобразуем активность в фактор
metadata_filtered$activity <- factor(metadata_filtered[[activity_col]])

# Оставляем только A0 и A1
metadata_filtered <- metadata_filtered[metadata_filtered$activity %in% c("A0", "A1"), ]

# Преобразуем в фактор (A0 - референс)
metadata_filtered$activity <- factor(metadata_filtered$activity, levels = c("A0", "A1"))

# basic cleaning
metadata_filtered$patient <- as.character(metadata_filtered$patient)
metadata_filtered$age <- as.numeric(as.character(metadata_filtered$age))


# ==============================================================================
# PART 3: RENAME DRUG VARIABLES
# ==============================================================================

colnames(metadata_filtered)[colnames(metadata_filtered) == "drug_1"] <- "5-ASA"
colnames(metadata_filtered)[colnames(metadata_filtered) == "drug_2"] <- "Glucocorticoids"
colnames(metadata_filtered)[colnames(metadata_filtered) == "drug_3"] <- "Immunosuppressors"
colnames(metadata_filtered)[colnames(metadata_filtered) == "drug_4"] <- "Genetically_Engineered_Therapies"
colnames(metadata_filtered)[colnames(metadata_filtered) == "drug_5"] <- "Inhibitor_JAK_kinase"
colnames(metadata_filtered)[colnames(metadata_filtered) == "drug_6"] <- "Antibiotics"

# ensure numeric
drug_cols <- c("5-ASA", "Glucocorticoids", "Immunosuppressors", 
               "Genetically_Engineered_Therapies", "Inhibitor_JAK_kinase", "Antibiotics")
metadata_filtered[drug_cols] <- lapply(metadata_filtered[drug_cols], function(x) as.numeric(as.character(x)))
metadata_filtered <- metadata_filtered[complete.cases(metadata_filtered[drug_cols]), ]
dim(metadata_filtered)

# ==============================================================================
# PART 4: ALIGN DATA
# ==============================================================================

common_samples <- intersect(rownames(metadata_filtered), rownames(input_data))

metadata_filtered <- metadata_filtered[common_samples, , drop = FALSE]
data_filtered <- input_data[common_samples, , drop = FALSE]

metadata_filtered$sex <- factor(metadata_filtered$sex, levels = c("F", "M"))

# Выбираем референсную группу (обычно A0 - неактивная)
activity_levels <- levels(metadata_filtered$activity)
if("A0" %in% activity_levels) {
  ref_level <- "A0"
} else {
  ref_level <- activity_levels[1]  # первый уровень как референс
}

cat("\nReference level for activity:", ref_level, "\n")

# ==============================================================================
# PART 5: RUN MAASLIN3 (ACTIVITY AS PREDICTOR)
# ==============================================================================

cat("\nRunning MaAsLin3 CD activity analysis...\n")

fit_data <- maaslin3(
  input_data     = data_filtered,
  input_metadata = metadata_filtered,
  output         = "Desktop/WORK/gut/1_stage/R/maaslin3_results_CD_activity1",
  
  fixed_effects  = c(
    "activity",
    "age",
    "sex",
    "5-ASA",
    "Glucocorticoids",
    "Immunosuppressors", 
    "Genetically_Engineered_Therapies",
    "Inhibitor_JAK_kinase",
    "Antibiotics"
  ),
  random_effects = c("patient"),
  
  reference = c(
    "activity=A0",
    "sex=F"
  ),
  
  normalization = "TSS",
  transform     = "LOG",
  correction    = "BH",
  
  min_abundance  = 0.0001,
  min_prevalence = 0.1,
  small_random_effects = TRUE
)


#==============================================================================
# PART 6: LOAD RESULTS
# ==============================================================================

results <- read.delim(
  "Desktop/WORK/gut/1_stage/R/maaslin3_results_CD_activity1/all_results.tsv",
  check.names = FALSE
)

# keep significant for er (main effect)
res_sig <- results %>% filter(pval_joint < 0.05)
res_cd_activity1 <- res_sig %>% filter(metadata == "activity")

# ==============================================================================
# PART 7: PLOT (ABUNDANCE) - Activity
# ==============================================================================

abund <- res_cd_activity1 %>%
  filter(model == "abundance") %>%
  filter(!is.na(coef), !is.na(stderr), !is.na(pval_joint)) %>%
  mutate(
    feature = fct_reorder(feature, coef),
    direction = ifelse(coef > 0, "Higher in A1", "Higher in A0"),
    lower = coef - 1.96 * stderr,
    upper = coef + 1.96 * stderr,
    neglogp = -log10(pval_joint)
  )

p_abund <- ggplot(abund, aes(x = coef, y = feature, color = direction)) +
  geom_point(aes(size = neglogp), alpha = 0.8) +
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.2, alpha = 0.6) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_color_manual(
    name = "Direction",
    values = c("Higher in A1" = "#1B4F72", "Higher in A0" = "#C06014")
  ) +
  scale_size_continuous(
    name = "-log10(P-value)"
  ) +
  scale_y_discrete(
    labels = function(x) gsub("_", " ", x)
  ) +
  labs(
    title = "Abundance associations",
    subtitle = "CD patients only: A1 vs A0\nadjusted for age, gender and therapy",
    x = "β coefficient (±95% CI)", 
    y = ""
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    legend.box = "vertical",
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    axis.text.y = element_text(face = "italic")
  )+
  guides(
    color = guide_legend(order = 1),
    size = guide_legend(order = 2)
  )

# ==============================================================================
# PART 8: PLOT (PREVALENCE)
# ==============================================================================

prev <- res_cd_activity1 %>%
  filter(model == "prevalence") %>%
  mutate(
    feature = fct_reorder(feature, coef),
    direction = ifelse(coef > 0, "More frequent in A1", "More frequent in A0"),
    lower = coef - 1.96 * stderr,
    upper = coef + 1.96 * stderr,
    neglogp = -log10(pval_joint)
  )

p_prev <- ggplot(prev, aes(x = coef, y = feature, color = direction)) +
  geom_point(aes(size = neglogp), alpha = 0.8) +
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.2, alpha = 0.6) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_color_manual(
    name = "Direction",
    values = c("More frequent in A1" = "#2E8B57", "More frequent in A0" = "#B22222")
  ) +
  scale_size_continuous(
    name = "-log10(P-value)"
  ) +
  scale_y_discrete(
    labels = function(x) gsub("_", " ", x)
  ) +
  labs(
    title = "Prevalence associations",
    subtitle = "CD patients only: A1 vs A0\nadjusted for age, gender and therapy",
    x = "β coefficient (±95% CI)", 
    y = ""
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    legend.box = "vertical",
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    axis.text.y = element_text(face = "italic")
  )+
  guides(
    color = guide_legend(order = 1),  # Direction ПЕРВАЯ
    size = guide_legend(order = 2)    # -log10(P-value) ВТОРАЯ
  )

# ==============================================================================
# PART 9: COMBINE & SAVE
# ==============================================================================

combined6 <- (p_abund | p_prev) +
  plot_annotation(
    title = "Bacterial biomarkers associated with A1 / A0 activity in Crohn's disease patients",
    theme = theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5, margin = margin(b = 15))
    )
  )







###################PART7




# ==============================================================================
# PART 2: FILTER CD PATIENTS AND CREATE ACTIVITY VARIABLE
# ==============================================================================

# Оставляем только CD
metadata_filtered <- input_metadata[input_metadata$diagnosis == "CD", ]

# Преобразуем активность в фактор
metadata_filtered$activity <- factor(metadata_filtered[[activity_col]])

# Оставляем только A0 и A2
metadata_filtered <- metadata_filtered[metadata_filtered$activity %in% c("A0", "A2"), ]

# Преобразуем в фактор (A0 - референс)
metadata_filtered$activity <- factor(metadata_filtered$activity, levels = c("A0", "A2"))

# basic cleaning
metadata_filtered$patient <- as.character(metadata_filtered$patient)
metadata_filtered$age <- as.numeric(as.character(metadata_filtered$age))


# ==============================================================================
# PART 3: RENAME DRUG VARIABLES
# ==============================================================================

colnames(metadata_filtered)[colnames(metadata_filtered) == "drug_1"] <- "5-ASA"
colnames(metadata_filtered)[colnames(metadata_filtered) == "drug_2"] <- "Glucocorticoids"
colnames(metadata_filtered)[colnames(metadata_filtered) == "drug_3"] <- "Immunosuppressors"
colnames(metadata_filtered)[colnames(metadata_filtered) == "drug_4"] <- "Genetically_Engineered_Therapies"
colnames(metadata_filtered)[colnames(metadata_filtered) == "drug_5"] <- "Inhibitor_JAK_kinase"
colnames(metadata_filtered)[colnames(metadata_filtered) == "drug_6"] <- "Antibiotics"

# ensure numeric
drug_cols <- c("5-ASA", "Glucocorticoids", "Immunosuppressors", 
               "Genetically_Engineered_Therapies", "Inhibitor_JAK_kinase", "Antibiotics")
metadata_filtered[drug_cols] <- lapply(metadata_filtered[drug_cols], function(x) as.numeric(as.character(x)))
metadata_filtered <- metadata_filtered[complete.cases(metadata_filtered[drug_cols]), ]
dim(metadata_filtered)

# ==============================================================================
# PART 4: ALIGN DATA
# ==============================================================================

common_samples <- intersect(rownames(metadata_filtered), rownames(input_data))

metadata_filtered <- metadata_filtered[common_samples, , drop = FALSE]
data_filtered <- input_data[common_samples, , drop = FALSE]

metadata_filtered$sex <- factor(metadata_filtered$sex, levels = c("F", "M"))

# Выбираем референсную группу (обычно A0 - неактивная)
activity_levels <- levels(metadata_filtered$activity)
if("A0" %in% activity_levels) {
  ref_level <- "A0"
} else {
  ref_level <- activity_levels[1]  # первый уровень как референс
}

cat("\nReference level for activity:", ref_level, "\n")

# ==============================================================================
# PART 5: RUN MAASLIN3 (ACTIVITY AS PREDICTOR)
# ==============================================================================

cat("\nRunning MaAsLin3 CD activity analysis...\n")

fit_data <- maaslin3(
  input_data     = data_filtered,
  input_metadata = metadata_filtered,
  output         = "Desktop/WORK/gut/1_stage/R/maaslin3_results_CD_activity2",
  
  fixed_effects  = c(
    "activity",
    "age",
    "sex",
    "5-ASA",
    "Glucocorticoids",
    "Immunosuppressors", 
    "Genetically_Engineered_Therapies",
    "Inhibitor_JAK_kinase",
    "Antibiotics"
  ),
  random_effects = c("patient"),
  
  reference = c(
    "activity=A0",
    "sex=F"
  ),
  
  normalization = "TSS",
  transform     = "LOG",
  correction    = "BH",
  
  min_abundance  = 0.0001,
  min_prevalence = 0.1,
  small_random_effects = TRUE
)


#==============================================================================
# PART 6: LOAD RESULTS
# ==============================================================================

results <- read.delim(
  "Desktop/WORK/gut/1_stage/R/maaslin3_results_CD_activity2/all_results.tsv",
  check.names = FALSE
)

# keep significant for er (main effect)
res_sig <- results %>% filter(pval_joint < 0.05)
res_cd_activity2 <- res_sig %>% filter(metadata == "activity")

# ==============================================================================
# PART 7: PLOT (ABUNDANCE) - Activity
# ==============================================================================

abund <- res_cd_activity2 %>%
  filter(model == "abundance") %>%
  filter(!is.na(coef), !is.na(stderr), !is.na(pval_joint)) %>%
  mutate(
    feature = fct_reorder(feature, coef),
    direction = ifelse(coef > 0, "Higher in A2", "Higher in A0"),
    lower = coef - 1.96 * stderr,
    upper = coef + 1.96 * stderr,
    neglogp = -log10(pval_joint)
  )

p_abund <- ggplot(abund, aes(x = coef, y = feature, color = direction)) +
  geom_point(aes(size = neglogp), alpha = 0.8) +
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.2, alpha = 0.6) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_color_manual(
    name = "Direction",
    values = c("Higher in A2" = "#1B4F72", "Higher in A0" = "#C06014")
  ) +
  scale_size_continuous(
    name = "-log10(P-value)"
  ) +
  scale_y_discrete(
    labels = function(x) gsub("_", " ", x)
  ) +
  labs(
    title = "Abundance associations",
    subtitle = "CD patients only: A2 vs A0\nadjusted for age, gender and therapy",
    x = "β coefficient (±95% CI)", 
    y = ""
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    legend.box = "vertical",
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    axis.text.y = element_text(face = "italic")
  )+
  guides(
    color = guide_legend(order = 1),  # Direction ПЕРВАЯ
    size = guide_legend(order = 2)    # -log10(P-value) ВТОРАЯ
  )

# ==============================================================================
# PART 8: PLOT (PREVALENCE)
# ==============================================================================

prev <- res_cd_activity2 %>%
  filter(model == "prevalence") %>%
  mutate(
    feature = fct_reorder(feature, coef),
    direction = ifelse(coef > 0, "More frequent in A2", "More frequent in A0"),
    lower = coef - 1.96 * stderr,
    upper = coef + 1.96 * stderr,
    neglogp = -log10(pval_joint)
  )

p_prev <- ggplot(prev, aes(x = coef, y = feature, color = direction)) +
  geom_point(aes(size = neglogp), alpha = 0.8) +
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.2, alpha = 0.6) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_color_manual(
    name = "Direction",
    values = c("More frequent in A2" = "#2E8B57", "More frequent in A0" = "#B22222")
  ) +
  scale_size_continuous(
    name = "-log10(P-value)"
  ) +
  scale_y_discrete(
    labels = function(x) gsub("_", " ", x)
  ) +
  labs(
    title = "Prevalence associations",
    subtitle = "CD patients only: A2 vs A0\nadjusted for age, gender and therapy",
    x = "β coefficient (±95% CI)", 
    y = ""
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    legend.box = "vertical",
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    axis.text.y = element_text(face = "italic")
  )+
  guides(
    color = guide_legend(order = 1),  # Direction ПЕРВАЯ
    size = guide_legend(order = 2)    # -log10(P-value) ВТОРАЯ
  )

# ==============================================================================
# PART 9: COMBINE & SAVE
# ==============================================================================

combined7 <- (p_abund | p_prev) +
  plot_annotation(
    title = "Bacterial biomarkers associated with A2 / A0 activity in Crohn's disease patients",
    theme = theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5, margin = margin(b = 15))
    )
  )




###################PART8




# ==============================================================================
# PART 2: FILTER UC PATIENTS AND CREATE ACTIVITY VARIABLE
# ==============================================================================

# Оставляем только UC
metadata_filtered <- input_metadata[input_metadata$diagnosis == "UC", ]

# Преобразуем активность в фактор
metadata_filtered$activity <- factor(metadata_filtered[[activity_col]])

# Оставляем только A0 и A1
metadata_filtered <- metadata_filtered[metadata_filtered$activity %in% c("A0", "A1"), ]

# Преобразуем в фактор (A0 - референс)
metadata_filtered$activity <- factor(metadata_filtered$activity, levels = c("A0", "A1"))

# basic cleaning
metadata_filtered$patient <- as.character(metadata_filtered$patient)
metadata_filtered$age <- as.numeric(as.character(metadata_filtered$age))


# ==============================================================================
# PART 3: RENAME DRUG VARIABLES
# ==============================================================================

colnames(metadata_filtered)[colnames(metadata_filtered) == "drug_1"] <- "5-ASA"
colnames(metadata_filtered)[colnames(metadata_filtered) == "drug_2"] <- "Glucocorticoids"
colnames(metadata_filtered)[colnames(metadata_filtered) == "drug_3"] <- "Immunosuppressors"
colnames(metadata_filtered)[colnames(metadata_filtered) == "drug_4"] <- "Genetically_Engineered_Therapies"
colnames(metadata_filtered)[colnames(metadata_filtered) == "drug_5"] <- "Inhibitor_JAK_kinase"
colnames(metadata_filtered)[colnames(metadata_filtered) == "drug_6"] <- "Antibiotics"

# ensure numeric
drug_cols <- c("5-ASA", "Glucocorticoids", "Immunosuppressors", 
               "Genetically_Engineered_Therapies", "Inhibitor_JAK_kinase", "Antibiotics")
metadata_filtered[drug_cols] <- lapply(metadata_filtered[drug_cols], function(x) as.numeric(as.character(x)))
metadata_filtered <- metadata_filtered[complete.cases(metadata_filtered[drug_cols]), ]
dim(metadata_filtered)

# ==============================================================================
# PART 4: ALIGN DATA
# ==============================================================================

common_samples <- intersect(rownames(metadata_filtered), rownames(input_data))

metadata_filtered <- metadata_filtered[common_samples, , drop = FALSE]
data_filtered <- input_data[common_samples, , drop = FALSE]

metadata_filtered$sex <- factor(metadata_filtered$sex, levels = c("F", "M"))

# Выбираем референсную группу (обычно A0 - неактивная)
activity_levels <- levels(metadata_filtered$activity)
if("A0" %in% activity_levels) {
  ref_level <- "A0"
} else {
  ref_level <- activity_levels[1]  # первый уровень как референс
}

cat("\nReference level for activity:", ref_level, "\n")

# ==============================================================================
# PART 5: RUN MAASLIN3 (ACTIVITY AS PREDICTOR)
# ==============================================================================

cat("\nRunning MaAsLin3 UC activity analysis...\n")

fit_data <- maaslin3(
  input_data     = data_filtered,
  input_metadata = metadata_filtered,
  output         = "Desktop/WORK/gut/1_stage/R/maaslin3_results_UC_activity1",
  
  fixed_effects  = c(
    "activity",
    "age",
    "sex",
    "5-ASA",
    "Glucocorticoids",
    "Immunosuppressors", 
    "Genetically_Engineered_Therapies",
    "Inhibitor_JAK_kinase"
  ),
  random_effects = c("patient"),
  
  reference = c(
    "activity=A0",
    "sex=F"
  ),
  
  normalization = "TSS",
  transform     = "LOG",
  correction    = "BH",
  
  min_abundance  = 0.0001,
  min_prevalence = 0.1,
  small_random_effects = TRUE
)


#==============================================================================
# PART 6: LOAD RESULTS
# ==============================================================================

results <- read.delim(
  "Desktop/WORK/gut/1_stage/R/maaslin3_results_UC_activity1/all_results.tsv",
  check.names = FALSE
)

# keep significant for er (main effect)
res_sig <- results %>% filter(pval_joint < 0.05)
res_uc_activity1 <- res_sig %>% filter(metadata == "activity")

# ==============================================================================
# PART 7: PLOT (ABUNDANCE) - Activity
# ==============================================================================

abund <- res_uc_activity1 %>%
  filter(model == "abundance") %>%
  filter(!is.na(coef), !is.na(stderr), !is.na(pval_joint)) %>%
  mutate(
    feature = fct_reorder(feature, coef),
    direction = ifelse(coef > 0, "Higher in A1", "Higher in A0"),
    lower = coef - 1.96 * stderr,
    upper = coef + 1.96 * stderr,
    neglogp = -log10(pval_joint)
  )

p_abund <- ggplot(abund, aes(x = coef, y = feature, color = direction)) +
  geom_point(aes(size = neglogp), alpha = 0.8) +
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.2, alpha = 0.6) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_color_manual(
    name = "Direction",
    values = c("Higher in A1" = "#1B4F72", "Higher in A0" = "#C06014")
  ) +
  scale_size_continuous(
    name = "-log10(P-value)"
  ) +
  scale_y_discrete(
    labels = function(x) gsub("_", " ", x)
  ) +
  labs(
    title = "Abundance associations",
    subtitle = "UC patients only: A1 vs A0\nadjusted for age, gender and therapy",
    x = "β coefficient (±95% CI)", 
    y = ""
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    legend.box = "vertical",
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    axis.text.y = element_text(face = "italic")
  )+
  guides(
    color = guide_legend(order = 1),
    size = guide_legend(order = 2)
  )

# ==============================================================================
# PART 8: PLOT (PREVALENCE)
# ==============================================================================

prev <- res_uc_activity1 %>%
  filter(model == "prevalence") %>%
  mutate(
    feature = fct_reorder(feature, coef),
    direction = ifelse(coef > 0, "More frequent in A1", "More frequent in A0"),
    lower = coef - 1.96 * stderr,
    upper = coef + 1.96 * stderr,
    neglogp = -log10(pval_joint)
  )

p_prev <- ggplot(prev, aes(x = coef, y = feature, color = direction)) +
  geom_point(aes(size = neglogp), alpha = 0.8) +
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.2, alpha = 0.6) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_color_manual(
    name = "Direction",
    values = c("More frequent in A1" = "#2E8B57", "More frequent in A0" = "#B22222")
  ) +
  scale_size_continuous(
    name = "-log10(P-value)"
  ) +
  scale_y_discrete(
    labels = function(x) gsub("_", " ", x)
  ) +
  labs(
    title = "Prevalence associations",
    subtitle = "UC patients only: A1 vs A0\nadjusted for age, gender and therapy",
    x = "β coefficient (±95% CI)", 
    y = ""
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    legend.box = "vertical",
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    axis.text.y = element_text(face = "italic")
  )+
  guides(
    color = guide_legend(order = 1),
    size = guide_legend(order = 2)
  )

# ==============================================================================
# PART 9: COMBINE & SAVE
# ==============================================================================

combined8 <- (p_abund | p_prev) +
  plot_annotation(
    title = "Bacterial biomarkers associated with A1 / A0 activity in Ulcerative Colitis patients",
    theme = theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5, margin = margin(b = 15))
    )
  )







###########################################################################PART9




# ==============================================================================
# PART 2: FILTER UC PATIENTS AND CREATE ACTIVITY VARIABLE
# ==============================================================================

# Оставляем только UC
metadata_filtered <- input_metadata[input_metadata$diagnosis == "UC", ]

# Преобразуем активность в фактор
metadata_filtered$activity <- factor(metadata_filtered[[activity_col]])

# Оставляем только A0 и A2
metadata_filtered <- metadata_filtered[metadata_filtered$activity %in% c("A0", "A2"), ]

# Преобразуем в фактор (A0 - референс)
metadata_filtered$activity <- factor(metadata_filtered$activity, levels = c("A0", "A2"))

# basic cleaning
metadata_filtered$patient <- as.character(metadata_filtered$patient)
metadata_filtered$age <- as.numeric(as.character(metadata_filtered$age))


# ==============================================================================
# PART 3: RENAME DRUG VARIABLES
# ==============================================================================

colnames(metadata_filtered)[colnames(metadata_filtered) == "drug_1"] <- "5-ASA"
colnames(metadata_filtered)[colnames(metadata_filtered) == "drug_2"] <- "Glucocorticoids"
colnames(metadata_filtered)[colnames(metadata_filtered) == "drug_3"] <- "Immunosuppressors"
colnames(metadata_filtered)[colnames(metadata_filtered) == "drug_4"] <- "Genetically_Engineered_Therapies"
colnames(metadata_filtered)[colnames(metadata_filtered) == "drug_5"] <- "Inhibitor_JAK_kinase"
colnames(metadata_filtered)[colnames(metadata_filtered) == "drug_6"] <- "Antibiotics"

# ensure numeric
drug_cols <- c("5-ASA", "Glucocorticoids", "Immunosuppressors", 
               "Genetically_Engineered_Therapies", "Inhibitor_JAK_kinase", "Antibiotics")
metadata_filtered[drug_cols] <- lapply(metadata_filtered[drug_cols], function(x) as.numeric(as.character(x)))
metadata_filtered <- metadata_filtered[complete.cases(metadata_filtered[drug_cols]), ]
dim(metadata_filtered)

# ==============================================================================
# PART 4: ALIGN DATA
# ==============================================================================

common_samples <- intersect(rownames(metadata_filtered), rownames(input_data))

metadata_filtered <- metadata_filtered[common_samples, , drop = FALSE]
data_filtered <- input_data[common_samples, , drop = FALSE]

metadata_filtered$sex <- factor(metadata_filtered$sex, levels = c("F", "M"))

# Выбираем референсную группу (обычно A0 - неактивная)
activity_levels <- levels(metadata_filtered$activity)
if("A0" %in% activity_levels) {
  ref_level <- "A0"
} else {
  ref_level <- activity_levels[1]  # первый уровень как референс
}

cat("\nReference level for activity:", ref_level, "\n")

# ==============================================================================
# PART 5: RUN MAASLIN3 (ACTIVITY AS PREDICTOR)
# ==============================================================================

cat("\nRunning MaAsLin3 UC activity analysis...\n")

fit_data <- maaslin3(
  input_data     = data_filtered,
  input_metadata = metadata_filtered,
  output         = "Desktop/WORK/gut/1_stage/R/maaslin3_results_UC_activity2",
  
  fixed_effects  = c(
    "activity",
    "age",
    "sex",
    "5-ASA",
    "Glucocorticoids",
    "Immunosuppressors", 
    "Genetically_Engineered_Therapies",
    "Inhibitor_JAK_kinase"
  ),
  random_effects = c("patient"),
  
  reference = c(
    "activity=A0",
    "sex=F"
  ),
  
  normalization = "TSS",
  transform     = "LOG",
  correction    = "BH",
  
  min_abundance  = 0.0001,
  min_prevalence = 0.1,
  small_random_effects = TRUE
)


#==============================================================================
# PART 6: LOAD RESULTS
# ==============================================================================

results <- read.delim(
  "Desktop/WORK/gut/1_stage/R/maaslin3_results_UC_activity2/all_results.tsv",
  check.names = FALSE
)

# keep significant for er (main effect)
res_sig <- results %>% filter(pval_joint < 0.05)
res_uc_activity2 <- res_sig %>% filter(metadata == "activity")

# ==============================================================================
# PART 7: PLOT (ABUNDANCE) - Activity
# ==============================================================================

abund <- res_uc_activity2 %>%
  filter(model == "abundance") %>%
  filter(!is.na(coef), !is.na(stderr), !is.na(pval_joint)) %>%
  mutate(
    feature = fct_reorder(feature, coef),
    direction = ifelse(coef > 0, "Higher in A2", "Higher in A0"),
    lower = coef - 1.96 * stderr,
    upper = coef + 1.96 * stderr,
    neglogp = -log10(pval_joint)
  )

p_abund <- ggplot(abund, aes(x = coef, y = feature, color = direction)) +
  geom_point(aes(size = neglogp), alpha = 0.8) +
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.2, alpha = 0.6) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_color_manual(
    name = "Direction",
    values = c("Higher in A2" = "#1B4F72", "Higher in A0" = "#C06014")
  ) +
  scale_size_continuous(
    name = "-log10(P-value)"
  ) +
  scale_y_discrete(
    labels = function(x) gsub("_", " ", x)
  ) +
  labs(
    title = "Abundance associations",
    subtitle = "UC patients only: A2 vs A0\nadjusted for age, gender and therapy",
    x = "β coefficient (±95% CI)", 
    y = ""
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    legend.box = "vertical",
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    axis.text.y = element_text(face = "italic")
  )+
  guides(
    color = guide_legend(order = 1),
    size = guide_legend(order = 2)
  )

# ==============================================================================
# PART 8: PLOT (PREVALENCE)
# ==============================================================================

prev <- res_uc_activity2 %>%
  filter(model == "prevalence") %>%
  mutate(
    feature = fct_reorder(feature, coef),
    direction = ifelse(coef > 0, "More frequent in A2", "More frequent in A0"),
    lower = coef - 1.96 * stderr,
    upper = coef + 1.96 * stderr,
    neglogp = -log10(pval_joint)
  )

p_prev <- ggplot(prev, aes(x = coef, y = feature, color = direction)) +
  geom_point(aes(size = neglogp), alpha = 0.8) +
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.2, alpha = 0.6) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_color_manual(
    name = "Direction",
    values = c("More frequent in A2" = "#2E8B57", "More frequent in A0" = "#B22222")
  ) +
  scale_size_continuous(
    name = "-log10(P-value)"
  ) +
  scale_y_discrete(
    labels = function(x) gsub("_", " ", x)
  ) +
  labs(
    title = "Prevalence associations",
    subtitle = "UC patients only: A2 vs A0\nadjusted for age, gender and therapy",
    x = "β coefficient (±95% CI)", 
    y = ""
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    legend.box = "vertical",
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    axis.text.y = element_text(face = "italic")
  )+
  guides(
    color = guide_legend(order = 1),
    size = guide_legend(order = 2)
  )

# ==============================================================================
# PART 9: COMBINE & SAVE
# ==============================================================================

combined9 <- (p_abund | p_prev) +
  plot_annotation(
    title = "Bacterial biomarkers associated with A2 / A0 activity in Ulcerative Colitis patients",
    theme = theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5, margin = margin(b = 15))
    )
  )



#####################################PART10



# ==============================================================================
# PART 2: FILTER CD PATIENTS AND CREATE Extraintestinal manifistation VARIABLE
# ==============================================================================

# Оставляем только CD
metadata_filtered <- input_metadata[input_metadata$diagnosis == "CD", ]



# ==============================================================================
# PART 3: RENAME DRUG VARIABLES
# ==============================================================================

colnames(metadata_filtered)[colnames(metadata_filtered) == "drug_1"] <- "5-ASA"
colnames(metadata_filtered)[colnames(metadata_filtered) == "drug_2"] <- "Glucocorticoids"
colnames(metadata_filtered)[colnames(metadata_filtered) == "drug_3"] <- "Immunosuppressors"
colnames(metadata_filtered)[colnames(metadata_filtered) == "drug_4"] <- "Genetically_Engineered_Therapies"
colnames(metadata_filtered)[colnames(metadata_filtered) == "drug_5"] <- "Inhibitor_JAK_kinase"
colnames(metadata_filtered)[colnames(metadata_filtered) == "drug_6"] <- "Antibiotics"

# ensure numeric
drug_cols <- c("5-ASA", "Glucocorticoids", "Immunosuppressors", 
               "Genetically_Engineered_Therapies", "Inhibitor_JAK_kinase", "Antibiotics")
metadata_filtered[drug_cols] <- lapply(metadata_filtered[drug_cols], function(x) as.numeric(as.character(x)))
metadata_filtered <- metadata_filtered[complete.cases(metadata_filtered[drug_cols]), ]
dim(metadata_filtered)
# ==============================================================================
# PART 4: ALIGN DATA
# ==============================================================================

common_samples <- intersect(rownames(metadata_filtered), rownames(input_data))

metadata_filtered <- metadata_filtered[common_samples, , drop = FALSE]
data_filtered <- input_data[common_samples, , drop = FALSE]

metadata_filtered$sex <- factor(metadata_filtered$sex, levels = c("F", "M"))
metadata_filtered$er <- factor(metadata_filtered$er, levels = c(0, 1))
metadata_filtered$em <- factor(metadata_filtered$em, levels = c(0, 1))
metadata_filtered$patient <- as.character(metadata_filtered$patient)
metadata_filtered$age <- as.numeric(as.character(metadata_filtered$age))

# ==============================================================================
# PART 5: RUN MAASLIN3 (Extraintestinal manifistation 1 vs 0 within CD, adjusting for ER)
# ==============================================================================

cat("Running MaAsLin3 CD-only Extraintestinal manifistation  analysis...\n")

fit_data <- maaslin3(
  input_data     = data_filtered,
  input_metadata = metadata_filtered,
  output         = "Desktop/WORK/gut/1_stage/R/maaslin3_results_CD_EIM",
  
  fixed_effects  = c(
    "em",        # основная переменная интереса
    "er",        # вычитание/поправка на er
    "age",
    "sex",
    "5-ASA",
    "Glucocorticoids",
    "Immunosuppressors",
    "Genetically_Engineered_Therapies",
    "Inhibitor_JAK_kinase",
    "Antibiotics"
  ),
  
  random_effects = c("patient"),
  
  reference = c(
    "em=0",      # em=0 как референс
    "er=0",      # er=0 как референс
    "sex=F"
  ),
  
  normalization = "TSS",
  transform     = "LOG",
  correction    = "BH",
  
  min_abundance  = 0.0001,
  min_prevalence = 0.1,
  
  small_random_effects = TRUE
)

#==============================================================================
# PART 6: LOAD RESULTS
# ==============================================================================

results <- read.delim(
  "Desktop/WORK/gut/1_stage/R/maaslin3_results_CD_EIM/all_results.tsv",
  check.names = FALSE
)

# keep significant for er (main effect)
res_sig <- results %>% filter(pval_joint < 0.05)
res_cd_em <- res_sig %>% filter(metadata == "em")

# ==============================================================================
# PART 7: PLOT (ABUNDANCE) - EIM 1 vs 0
# ==============================================================================

abund <- res_cd_em %>%
  filter(model == "abundance") %>%
  filter(!is.na(coef), !is.na(stderr), !is.na(pval_joint)) %>%
  mutate(
    feature = fct_reorder(feature, coef),
    direction = ifelse(coef > 0, "Higher in eim (yes)", "Higher in eim (no)"),
    lower = coef - 1.96 * stderr,
    upper = coef + 1.96 * stderr,
    neglogp = -log10(pval_joint)
  )

p_abund <- ggplot(abund, aes(x = coef, y = feature, color = direction)) +
  geom_point(aes(size = neglogp), alpha = 0.8) +
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.2, alpha = 0.6) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_color_manual(
    name = "Direction",
    values = c("Higher in eim (yes)" = "#1B4F72", "Higher in eim (no)" = "#C06014")
  ) +
  scale_size_continuous(
    name = "-log10(P-value)"
  ) +
  scale_y_discrete(
    labels = function(x) gsub("_", " ", x)
  ) +
  labs(
    title = "Abundance associations",
    subtitle = "CD patients only, extraintestinal manifistation (eim): yes vs no \nadjusted for age, gender, therapy and exacerbation / remission",
    x = "β coefficient (±95% CI)", 
    y = ""
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    legend.box = "vertical",
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    axis.text.y = element_text(face = "italic")
  )+
  guides(
    color = guide_legend(order = 1),  # Direction ПЕРВАЯ
    size = guide_legend(order = 2)    # -log10(P-value) ВТОРАЯ
  )

# ==============================================================================
# PART 8: PLOT (PREVALENCE)
# ==============================================================================

prev <- res_cd_em %>%
  filter(model == "prevalence") %>%
  mutate(
    feature = fct_reorder(feature, coef),
    direction = ifelse(coef > 0, "More frequent in eim (yes)", "More frequent in eim (no)"),
    lower = coef - 1.96 * stderr,
    upper = coef + 1.96 * stderr,
    neglogp = -log10(pval_joint)
  )

p_prev <- ggplot(prev, aes(x = coef, y = feature, color = direction)) +
  geom_point(aes(size = neglogp), alpha = 0.8) +
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.2, alpha = 0.6) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_color_manual(
    name = "Direction",
    values = c("More frequent in eim (yes)" = "#2E8B57", "More frequent in eim (no)" = "#B22222")
  ) +
  scale_size_continuous(
    name = "-log10(P-value)"
  ) +
  scale_y_discrete(
    labels = function(x) gsub("_", " ", x)
  ) +
  labs(
    title = "Prevalence associations",
    subtitle = "CD patients only, extraintestinal manifistation (eim): yes vx no \nadjusted for age, gender, therapy and exacerbation / remission",
    x = "β coefficient (±95% CI)", 
    y = ""
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    legend.box = "vertical",
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    axis.text.y = element_text(face = "italic")
  )+
  guides(
    color = guide_legend(order = 1),  # Direction ПЕРВАЯ
    size = guide_legend(order = 2)    # -log10(P-value) ВТОРАЯ
  )

# ==============================================================================
# PART 9: COMBINE & SAVE
# ==============================================================================

combined10 <- (p_abund | p_prev) +
  plot_annotation(
    title = "Bacterial biomarkers associated with extraintestinal manifistation in Crohn's disease patients",
    theme = theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5, margin = margin(b = 15))
    )
  )



#####################################PART11



# ==============================================================================
# PART 2: FILTER CD PATIENTS AND CREATE Extraintestinal manifistation VARIABLE
# ==============================================================================

# Оставляем только CD
metadata_filtered <- input_metadata[input_metadata$diagnosis == "UC", ]


# ==============================================================================
# PART 3: RENAME DRUG VARIABLES
# ==============================================================================

colnames(metadata_filtered)[colnames(metadata_filtered) == "drug_1"] <- "5-ASA"
colnames(metadata_filtered)[colnames(metadata_filtered) == "drug_2"] <- "Glucocorticoids"
colnames(metadata_filtered)[colnames(metadata_filtered) == "drug_3"] <- "Immunosuppressors"
colnames(metadata_filtered)[colnames(metadata_filtered) == "drug_4"] <- "Genetically_Engineered_Therapies"
colnames(metadata_filtered)[colnames(metadata_filtered) == "drug_5"] <- "Inhibitor_JAK_kinase"
colnames(metadata_filtered)[colnames(metadata_filtered) == "drug_6"] <- "Antibiotics"

# ensure numeric
drug_cols <- c("5-ASA", "Glucocorticoids", "Immunosuppressors", 
               "Genetically_Engineered_Therapies", "Inhibitor_JAK_kinase", "Antibiotics")
metadata_filtered[drug_cols] <- lapply(metadata_filtered[drug_cols], function(x) as.numeric(as.character(x)))
metadata_filtered <- metadata_filtered[complete.cases(metadata_filtered[drug_cols]), ]
dim(metadata_filtered)
# ==============================================================================
# PART 4: ALIGN DATA
# ==============================================================================

common_samples <- intersect(rownames(metadata_filtered), rownames(input_data))

metadata_filtered <- metadata_filtered[common_samples, , drop = FALSE]
data_filtered <- input_data[common_samples, , drop = FALSE]

metadata_filtered$sex <- factor(metadata_filtered$sex, levels = c("F", "M"))
metadata_filtered$er <- factor(metadata_filtered$er, levels = c(0, 1))
metadata_filtered$em <- factor(metadata_filtered$em, levels = c(0, 1))
metadata_filtered$patient <- as.character(metadata_filtered$patient)
metadata_filtered$age <- as.numeric(as.character(metadata_filtered$age))

# ==============================================================================
# PART 5: RUN MAASLIN3 (Extraintestinal manifistation 1 vs 0 within CD, adjusting for ER)
# ==============================================================================

cat("Running MaAsLin3 CD-only Extraintestinal manifistation  analysis...\n")

fit_data <- maaslin3(
  input_data     = data_filtered,
  input_metadata = metadata_filtered,
  output         = "Desktop/WORK/gut/1_stage/R/maaslin3_results_UC_EIM",
  
  fixed_effects  = c(
    "em",        # основная переменная интереса
    "er",        # вычитание/поправка на er
    "age",
    "sex",
    "5-ASA",
    "Glucocorticoids",
    "Immunosuppressors",
    "Genetically_Engineered_Therapies",
    "Inhibitor_JAK_kinase"
  ),
  
  random_effects = c("patient"),
  
  reference = c(
    "em=0",      # em=0 как референс
    "er=0",      # er=0 как референс
    "sex=F"
  ),
  
  normalization = "TSS",
  transform     = "LOG",
  correction    = "BH",
  
  min_abundance  = 0.0001,
  min_prevalence = 0.1,
  
  small_random_effects = TRUE
)

#==============================================================================
# PART 6: LOAD RESULTS
# ==============================================================================

results <- read.delim(
  "Desktop/WORK/gut/1_stage/R/maaslin3_results_UC_EIM/all_results.tsv",
  check.names = FALSE
)

# keep significant for er (main effect)
res_sig <- results %>% filter(pval_joint < 0.05)
res_uc_em <- res_sig %>% filter(metadata == "em")

# ==============================================================================
# PART 7: PLOT (ABUNDANCE) - EIM 1 vs 0
# ==============================================================================

abund <- res_uc_em %>%
  filter(model == "abundance") %>%
  filter(!is.na(coef), !is.na(stderr), !is.na(pval_joint)) %>%
  mutate(
    feature = fct_reorder(feature, coef),
    direction = ifelse(coef > 0, "Higher in eim (yes)", "Higher in eim (no)"),
    lower = coef - 1.96 * stderr,
    upper = coef + 1.96 * stderr,
    neglogp = -log10(pval_joint)
  )

p_abund <- ggplot(abund, aes(x = coef, y = feature, color = direction)) +
  geom_point(aes(size = neglogp), alpha = 0.8) +
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.2, alpha = 0.6) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_color_manual(
    name = "Direction",
    values = c("Higher in eim (yes)" = "#1B4F72", "Higher in eim (no)" = "#C06014")
  ) +
  scale_size_continuous(
    name = "-log10(P-value)"
  ) +
  scale_y_discrete(
    labels = function(x) gsub("_", " ", x)
  ) +
  labs(
    title = "Abundance associations",
    subtitle = "UC patients only, extraintestinal manifistation (eim): yes vs no \nadjusted for age, gender, therapy and exacerbation / remission",
    x = "β coefficient (±95% CI)", 
    y = ""
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    legend.box = "vertical",
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    axis.text.y = element_text(face = "italic")
  )+
  guides(
    color = guide_legend(order = 1),  # Direction ПЕРВАЯ
    size = guide_legend(order = 2)    # -log10(P-value) ВТОРАЯ
  )

# ==============================================================================
# PART 8: PLOT (PREVALENCE)
# ==============================================================================

prev <- res_uc_em %>%
  filter(model == "prevalence") %>%
  mutate(
    feature = fct_reorder(feature, coef),
    direction = ifelse(coef > 0, "More frequent in eim (yes)", "More frequent in eim (no)"),
    lower = coef - 1.96 * stderr,
    upper = coef + 1.96 * stderr,
    neglogp = -log10(pval_joint)
  )

p_prev <- ggplot(prev, aes(x = coef, y = feature, color = direction)) +
  geom_point(aes(size = neglogp), alpha = 0.8) +
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.2, alpha = 0.6) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_color_manual(
    name = "Direction",
    values = c("More frequent in eim (yes)" = "#2E8B57", "More frequent in eim (no)" = "#B22222")
  ) +
  scale_size_continuous(
    name = "-log10(P-value)"
  ) +
  scale_y_discrete(
    labels = function(x) gsub("_", " ", x)
  ) +
  labs(
    title = "Prevalence associations",
    subtitle = "UC patients only, extraintestinal manifistation (eim): yes vx no \nadjusted for age, gender, therapy and exacerbation / remission",
    x = "β coefficient (±95% CI)", 
    y = ""
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    legend.box = "vertical",
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    axis.text.y = element_text(face = "italic")
  )+
  guides(
    color = guide_legend(order = 1),  # Direction ПЕРВАЯ
    size = guide_legend(order = 2)    # -log10(P-value) ВТОРАЯ
  )

# ==============================================================================
# PART 9: COMBINE & SAVE
# ==============================================================================

combined11 <- (p_abund | p_prev) +
  plot_annotation(
    title = "Bacterial biomarkers associated with extraintestinal manifistation in Ulcerative Colitis patients",
    theme = theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5, margin = margin(b = 15))
    )
  )



###################PART12




# ==============================================================================
# PART 2: FILTER CD PATIENTS AND CREATE EF VARIABLE
# ==============================================================================

# Оставляем только CD
metadata_filtered <- input_metadata[input_metadata$diagnosis == "CD", ]

# Оставляем только L1 и L2
metadata_filtered <- metadata_filtered[metadata_filtered$ef %in% c("L1", "L2"), ]

# Преобразуем в фактор (L1 - референс)
metadata_filtered$ef <- factor(metadata_filtered$ef, levels = c("L1", "L2"))



# ==============================================================================
# PART 3: RENAME DRUG VARIABLES
# ==============================================================================

colnames(metadata_filtered)[colnames(metadata_filtered) == "drug_1"] <- "5-ASA"
colnames(metadata_filtered)[colnames(metadata_filtered) == "drug_2"] <- "Glucocorticoids"
colnames(metadata_filtered)[colnames(metadata_filtered) == "drug_3"] <- "Immunosuppressors"
colnames(metadata_filtered)[colnames(metadata_filtered) == "drug_4"] <- "Genetically_Engineered_Therapies"
colnames(metadata_filtered)[colnames(metadata_filtered) == "drug_5"] <- "Inhibitor_JAK_kinase"
colnames(metadata_filtered)[colnames(metadata_filtered) == "drug_6"] <- "Antibiotics"

# ensure numeric
drug_cols <- c("5-ASA", "Glucocorticoids", "Immunosuppressors", 
               "Genetically_Engineered_Therapies", "Inhibitor_JAK_kinase", "Antibiotics")
metadata_filtered[drug_cols] <- lapply(metadata_filtered[drug_cols], function(x) as.numeric(as.character(x)))
metadata_filtered <- metadata_filtered[complete.cases(metadata_filtered[drug_cols]), ]
dim(metadata_filtered)

# ==============================================================================
# PART 4: ALIGN DATA
# ==============================================================================

common_samples <- intersect(rownames(metadata_filtered), rownames(input_data))

metadata_filtered <- metadata_filtered[common_samples, , drop = FALSE]
data_filtered <- input_data[common_samples, , drop = FALSE]

metadata_filtered$sex <- factor(metadata_filtered$sex, levels = c("F", "M"))
metadata_filtered$patient <- as.character(metadata_filtered$patient)
metadata_filtered$age <- as.numeric(as.character(metadata_filtered$age))
metadata_filtered$er <- factor(metadata_filtered$er, levels = c(0, 1))
metadata_filtered$em <- factor(metadata_filtered$em, levels = c(0, 1))

# ==============================================================================
# PART 5: RUN MAASLIN3 (ACTIVITY AS PREDICTOR)
# ==============================================================================

cat("\nRunning MaAsLin3 CD activity analysis...\n")

fit_data <- maaslin3(
  input_data     = data_filtered,
  input_metadata = metadata_filtered,
  output         = "Desktop/WORK/gut/1_stage/R/maaslin3_results_CD_ef1",
  
  fixed_effects  = c(
    "ef",
    "er",
    "em",
    "age",
    "sex",
    "5-ASA",
    "Glucocorticoids",
    "Immunosuppressors", 
    "Genetically_Engineered_Therapies",
    "Inhibitor_JAK_kinase"
  ),
  
  random_effects = c("patient"),
  
  reference = c(
    "ef=L1",
    "er=0",
    "em=0",
    "sex=F"
  ),
  
  normalization = "TSS",
  transform     = "LOG",
  correction    = "BH",
  
  min_abundance  = 0.0001,
  min_prevalence = 0.1,
  small_random_effects = TRUE
)

#==============================================================================
# PART 6: LOAD RESULTS
# ==============================================================================

results <- read.delim(
  "Desktop/WORK/gut/1_stage/R/maaslin3_results_CD_ef1/all_results.tsv",
  check.names = FALSE
)

# keep significant for er (main effect)
res_sig <- results %>% filter(pval_joint < 0.05)
res_cd_ef1 <- res_sig %>% filter(metadata == "ef")

# ==============================================================================
# PART 7: PLOT (ABUNDANCE) - Activity
# ==============================================================================

abund <- res_cd_ef1 %>%
  filter(model == "abundance") %>%
  filter(!is.na(coef), !is.na(stderr), !is.na(pval_joint)) %>%
  mutate(
    feature = fct_reorder(feature, coef),
    direction = ifelse(coef > 0, "Higher in L2", "Higher in L1"),
    lower = coef - 1.96 * stderr,
    upper = coef + 1.96 * stderr,
    neglogp = -log10(pval_joint)
  )

p_abund <- ggplot(abund, aes(x = coef, y = feature, color = direction)) +
  geom_point(aes(size = neglogp), alpha = 0.8) +
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.2, alpha = 0.6) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_color_manual(
    name = "Direction",
    values = c("Higher in L2" = "#1B4F72", "Higher in L1" = "#C06014")
  ) +
  scale_size_continuous(
    name = "-log10(P-value)"
  ) +
  scale_y_discrete(
    labels = function(x) gsub("_", " ", x)
  ) +
  labs(
    title = "Abundance associations",
    subtitle = "CD patients only, endoscopic findings: L2 vs L1\nadjusted for age, gender, therapy, exacerbation / remission and\nextraintestinal manifistation",
    x = "β coefficient (±95% CI)", 
    y = ""
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    legend.box = "vertical",
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    axis.text.y = element_text(face = "italic")
  )+
  guides(
    color = guide_legend(order = 1),  # Direction ПЕРВАЯ
    size = guide_legend(order = 2)    # -log10(P-value) ВТОРАЯ
  )

# ==============================================================================
# PART 8: PLOT (PREVALENCE)
# ==============================================================================

prev <- res_cd_ef1 %>%
  filter(model == "prevalence") %>%
  mutate(
    feature = fct_reorder(feature, coef),
    direction = ifelse(coef > 0, "More frequent in L2", "More frequent in L1"),
    lower = coef - 1.96 * stderr,
    upper = coef + 1.96 * stderr,
    neglogp = -log10(pval_joint)
  )

p_prev <- ggplot(prev, aes(x = coef, y = feature, color = direction)) +
  geom_point(aes(size = neglogp), alpha = 0.8) +
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.2, alpha = 0.6) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_color_manual(
    name = "Direction",
    values = c("More frequent in L2" = "#2E8B57", "More frequent in L1" = "#B22222")
  ) +
  scale_size_continuous(
    name = "-log10(P-value)"
  ) +
  scale_y_discrete(
    labels = function(x) gsub("_", " ", x)
  ) +
  labs(
    title = "Prevalence associations",
    subtitle = "CD patients only, endoscopic findings: L2 vs L1\nadjusted for age, gender, therapy, exacerbation / remission and\nextraintestinal manifistation",
    x = "β coefficient (±95% CI)", 
    y = ""
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    legend.box = "vertical",
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    axis.text.y = element_text(face = "italic")
  )+
  guides(
    color = guide_legend(order = 1),  # Direction ПЕРВАЯ
    size = guide_legend(order = 2)    # -log10(P-value) ВТОРАЯ
  )

# ==============================================================================
# PART 9: COMBINE & SAVE
# ==============================================================================

combined12 <- (p_abund | p_prev) +
  plot_annotation(
    title = "Bacterial biomarkers associated with endoscopic findings (L2 / L1) in Crohn's disease patients",
    theme = theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5, margin = margin(b = 15))
    )
  )



###################PART13




# ==============================================================================
# PART 2: FILTER CD PATIENTS AND CREATE EF VARIABLE
# ==============================================================================

# Оставляем только CD
metadata_filtered <- input_metadata[input_metadata$diagnosis == "CD", ]

# Оставляем только L1 и L3
metadata_filtered <- metadata_filtered[metadata_filtered$ef %in% c("L1", "L3"), ]

# Преобразуем в фактор (L1 - референс)
metadata_filtered$ef<- factor(metadata_filtered$ef, levels = c("L1", "L3"))

# ==============================================================================
# PART 3: RENAME DRUG VARIABLES
# ==============================================================================

colnames(metadata_filtered)[colnames(metadata_filtered) == "drug_1"] <- "5-ASA"
colnames(metadata_filtered)[colnames(metadata_filtered) == "drug_2"] <- "Glucocorticoids"
colnames(metadata_filtered)[colnames(metadata_filtered) == "drug_3"] <- "Immunosuppressors"
colnames(metadata_filtered)[colnames(metadata_filtered) == "drug_4"] <- "Genetically_Engineered_Therapies"
colnames(metadata_filtered)[colnames(metadata_filtered) == "drug_5"] <- "Inhibitor_JAK_kinase"
colnames(metadata_filtered)[colnames(metadata_filtered) == "drug_6"] <- "Antibiotics"

# ensure numeric
drug_cols <- c("5-ASA", "Glucocorticoids", "Immunosuppressors", 
               "Genetically_Engineered_Therapies", "Inhibitor_JAK_kinase", "Antibiotics")
metadata_filtered[drug_cols] <- lapply(metadata_filtered[drug_cols], function(x) as.numeric(as.character(x)))
metadata_filtered <- metadata_filtered[complete.cases(metadata_filtered[drug_cols]), ]
dim(metadata_filtered)

# ==============================================================================
# PART 4: ALIGN DATA
# ==============================================================================

common_samples <- intersect(rownames(metadata_filtered), rownames(input_data))

metadata_filtered <- metadata_filtered[common_samples, , drop = FALSE]
data_filtered <- input_data[common_samples, , drop = FALSE]

metadata_filtered$sex <- factor(metadata_filtered$sex, levels = c("F", "M"))
metadata_filtered$patient <- as.character(metadata_filtered$patient)
metadata_filtered$age <- as.numeric(as.character(metadata_filtered$age))
metadata_filtered$er <- factor(metadata_filtered$er, levels = c(0, 1))
metadata_filtered$em <- factor(metadata_filtered$em, levels = c(0, 1))

# ==============================================================================
# PART 5: RUN MAASLIN3 (ACTIVITY AS PREDICTOR)
# ==============================================================================

cat("\nRunning MaAsLin3 CD activity analysis...\n")

fit_data <- maaslin3(
  input_data     = data_filtered,
  input_metadata = metadata_filtered,
  output         = "Desktop/WORK/gut/1_stage/R/maaslin3_results_CD_ef2",
  
  fixed_effects  = c(
    "ef",
    "er",
    "em",
    "age",
    "sex",
    "5-ASA",
    "Glucocorticoids",
    "Immunosuppressors", 
    "Genetically_Engineered_Therapies",
    "Inhibitor_JAK_kinase",
    "Antibiotics"
  ),
  
  random_effects = c("patient"),
  
  reference = c(
    "ef=L1",
    "er=0",
    "em=0",
    "sex=F"
  ),
  
  normalization = "TSS",
  transform     = "LOG",
  correction    = "BH",
  
  min_abundance  = 0.0001,
  min_prevalence = 0.1,
  small_random_effects = TRUE
)


#==============================================================================
# PART 6: LOAD RESULTS
# ==============================================================================

results <- read.delim(
  "Desktop/WORK/gut/1_stage/R/maaslin3_results_CD_ef2/all_results.tsv",
  check.names = FALSE
)

# keep significant for er (main effect)
res_sig <- results %>% filter(pval_joint < 0.05)
res_cd_ef2 <- res_sig %>% filter(metadata == "ef")

# ==============================================================================
# PART 7: PLOT (ABUNDANCE) - Activity
# ==============================================================================

abund <- res_cd_ef2 %>%
  filter(model == "abundance") %>%
  filter(!is.na(coef), !is.na(stderr), !is.na(pval_joint)) %>%
  mutate(
    feature = fct_reorder(feature, coef),
    direction = ifelse(coef > 0, "Higher in L3", "Higher in L1"),
    lower = coef - 1.96 * stderr,
    upper = coef + 1.96 * stderr,
    neglogp = -log10(pval_joint)
  )

p_abund <- ggplot(abund, aes(x = coef, y = feature, color = direction)) +
  geom_point(aes(size = neglogp), alpha = 0.8) +
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.2, alpha = 0.6) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_color_manual(
    name = "Direction",
    values = c("Higher in L3" = "#1B4F72", "Higher in L1" = "#C06014")
  ) +
  scale_size_continuous(
    name = "-log10(P-value)"
  ) +
  scale_y_discrete(
    labels = function(x) gsub("_", " ", x)
  ) +
  labs(
    title = "Abundance associations",
    subtitle = "CD patients only, endoscopic findings: L3 vs L1\nadjusted for age, gender, therapy, exacerbation / remission and\nextraintestinal manifistation",
    x = "β coefficient (±95% CI)", 
    y = ""
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    legend.box = "vertical",
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    axis.text.y = element_text(face = "italic")
  )+
  guides(
    color = guide_legend(order = 1),  # Direction ПЕРВАЯ
    size = guide_legend(order = 2)    # -log10(P-value) ВТОРАЯ
  )

# ==============================================================================
# PART 8: PLOT (PREVALENCE)
# ==============================================================================

prev <- res_cd_ef2 %>%
  filter(model == "prevalence") %>%
  mutate(
    feature = fct_reorder(feature, coef),
    direction = ifelse(coef > 0, "More frequent in L3", "More frequent in L1"),
    lower = coef - 1.96 * stderr,
    upper = coef + 1.96 * stderr,
    neglogp = -log10(pval_joint)
  )

p_prev <- ggplot(prev, aes(x = coef, y = feature, color = direction)) +
  geom_point(aes(size = neglogp), alpha = 0.8) +
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.2, alpha = 0.6) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_color_manual(
    name = "Direction",
    values = c("More frequent in L3" = "#2E8B57", "More frequent in L1" = "#B22222")
  ) +
  scale_size_continuous(
    name = "-log10(P-value)"
  ) +
  scale_y_discrete(
    labels = function(x) gsub("_", " ", x)
  ) +
  labs(
    title = "Prevalence associations",
    subtitle = "CD patients only, endoscopic findings: L3 vs L1\nadjusted for age, gender, therapy, exacerbation / remission and\nextraintestinal manifistation",
    x = "β coefficient (±95% CI)", 
    y = ""
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    legend.box = "vertical",
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    axis.text.y = element_text(face = "italic")
  )+
  guides(
    color = guide_legend(order = 1),  # Direction ПЕРВАЯ
    size = guide_legend(order = 2)    # -log10(P-value) ВТОРАЯ
  )

# ==============================================================================
# PART 9: COMBINE & SAVE
# ==============================================================================

combined13 <- (p_abund | p_prev) +
  plot_annotation(
    title = "Bacterial biomarkers associated with endoscopic findings (L3 / L1) in Crohn's disease patients",
    theme = theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5, margin = margin(b = 15))
    )
  )






##########################################################################PART14




# ==============================================================================
# PART 2: FILTER UC PATIENTS AND CREATE EF VARIABLE
# ==============================================================================

# Оставляем только UC
metadata_filtered <- input_metadata[input_metadata$diagnosis == "UC", ]

# Оставляем только E2 и E3
metadata_filtered <- metadata_filtered[metadata_filtered$ef %in% c("E2", "E3"), ]

# Преобразуем в фактор (E2 - референс)
metadata_filtered$ef <- factor(metadata_filtered$ef, levels = c("E2", "E3"))

# ==============================================================================
# PART 3: RENAME DRUG VARIABLES
# ==============================================================================

colnames(metadata_filtered)[colnames(metadata_filtered) == "drug_1"] <- "5-ASA"
colnames(metadata_filtered)[colnames(metadata_filtered) == "drug_2"] <- "Glucocorticoids"
colnames(metadata_filtered)[colnames(metadata_filtered) == "drug_3"] <- "Immunosuppressors"
colnames(metadata_filtered)[colnames(metadata_filtered) == "drug_4"] <- "Genetically_Engineered_Therapies"
colnames(metadata_filtered)[colnames(metadata_filtered) == "drug_5"] <- "Inhibitor_JAK_kinase"
colnames(metadata_filtered)[colnames(metadata_filtered) == "drug_6"] <- "Antibiotics"

# ensure numeric
drug_cols <- c("5-ASA", "Glucocorticoids", "Immunosuppressors", 
               "Genetically_Engineered_Therapies", "Inhibitor_JAK_kinase", "Antibiotics")
metadata_filtered[drug_cols] <- lapply(metadata_filtered[drug_cols], function(x) as.numeric(as.character(x)))
metadata_filtered <- metadata_filtered[complete.cases(metadata_filtered[drug_cols]), ]
dim(metadata_filtered)

# ==============================================================================
# PART 4: ALIGN DATA
# ==============================================================================

common_samples <- intersect(rownames(metadata_filtered), rownames(input_data))

metadata_filtered <- metadata_filtered[common_samples, , drop = FALSE]
data_filtered <- input_data[common_samples, , drop = FALSE]

metadata_filtered$sex <- factor(metadata_filtered$sex, levels = c("F", "M"))
metadata_filtered$patient <- as.character(metadata_filtered$patient)
metadata_filtered$age <- as.numeric(as.character(metadata_filtered$age))
metadata_filtered$er <- factor(metadata_filtered$er, levels = c(0, 1))
metadata_filtered$em <- factor(metadata_filtered$em, levels = c(0, 1))
# ==============================================================================
# PART 5: RUN MAASLIN3 (ACTIVITY AS PREDICTOR)
# ==============================================================================

cat("\nRunning MaAsLin3 UC activity analysis...\n")

fit_data <- maaslin3(
  input_data     = data_filtered,
  input_metadata = metadata_filtered,
  output         = "Desktop/WORK/gut/1_stage/R/maaslin3_results_UC_ef1",
  
  fixed_effects  = c(
    "ef",
    "er",
    "em",
    "age",
    "sex",
    "5-ASA",
    "Glucocorticoids",
    "Immunosuppressors", 
    "Genetically_Engineered_Therapies",
    "Inhibitor_JAK_kinase"
  ),
  
  random_effects = c("patient"),
  
  reference = c(
    "ef=E2",
    "er=0",
    "em=0",
    "sex=F"
  ),
  
  normalization = "TSS",
  transform     = "LOG",
  correction    = "BH",
  
  min_abundance  = 0.0001,
  min_prevalence = 0.1,
  small_random_effects = TRUE
)


#==============================================================================
# PART 6: LOAD RESULTS
# ==============================================================================

results <- read.delim(
  "Desktop/WORK/gut/1_stage/R/maaslin3_results_UC_ef1/all_results.tsv",
  check.names = FALSE
)

# keep significant for er (main effect)
res_sig <- results %>% filter(pval_joint < 0.05)
res_uc_ef1 <- res_sig %>% filter(metadata == "ef")

# ==============================================================================
# PART 7: PLOT (ABUNDANCE) - Activity
# ==============================================================================

abund <- res_uc_ef1 %>%
  filter(model == "abundance") %>%
  filter(!is.na(coef), !is.na(stderr), !is.na(pval_joint)) %>%
  mutate(
    feature = fct_reorder(feature, coef),
    direction = ifelse(coef > 0, "Higher in E3", "Higher in E2"),
    lower = coef - 1.96 * stderr,
    upper = coef + 1.96 * stderr,
    neglogp = -log10(pval_joint)
  )

p_abund <- ggplot(abund, aes(x = coef, y = feature, color = direction)) +
  geom_point(aes(size = neglogp), alpha = 0.8) +
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.2, alpha = 0.6) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_color_manual(
    name = "Direction",
    values = c("Higher in E3" = "#1B4F72", "Higher in E2" = "#C06014")
  ) +
  scale_size_continuous(
    name = "-log10(P-value)"
  ) +
  scale_y_discrete(
    labels = function(x) gsub("_", " ", x)
  ) +
  labs(
    title = "Abundance associations",
    subtitle = "UC patients only, endoscopic findings: E3 vs E2\nadjusted for age, gender, therapy, exacerbation / remission and\nextraintestinal manifistation",
    x = "β coefficient (±95% CI)", 
    y = ""
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    legend.box = "vertical",
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    axis.text.y = element_text(face = "italic")
  )+
  guides(
    color = guide_legend(order = 1),  # Direction ПЕРВАЯ
    size = guide_legend(order = 2)    # -log10(P-value) ВТОРАЯ
  )

# ==============================================================================
# PART 8: PLOT (PREVALENCE)
# ==============================================================================

prev <- res_uc_ef1 %>%
  filter(model == "prevalence") %>%
  mutate(
    feature = fct_reorder(feature, coef),
    direction = ifelse(coef > 0, "More frequent in E3", "More frequent in E2"),
    lower = coef - 1.96 * stderr,
    upper = coef + 1.96 * stderr,
    neglogp = -log10(pval_joint)
  )

p_prev <- ggplot(prev, aes(x = coef, y = feature, color = direction)) +
  geom_point(aes(size = neglogp), alpha = 0.8) +
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.2, alpha = 0.6) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_color_manual(
    name = "Direction",
    values = c("More frequent in E3" = "#2E8B57", "More frequent in E2" = "#B22222")
  ) +
  scale_size_continuous(
    name = "-log10(P-value)"
  ) +
  scale_y_discrete(
    labels = function(x) gsub("_", " ", x)
  ) +
  labs(
    title = "Prevalence associations",
    subtitle = "UC patients only, endoscopic findings: E3 vs E2\nadjusted for age, gender, therapy, exacerbation / remission and\nextraintestinal manifistation",
    x = "β coefficient (±95% CI)", 
    y = ""
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    legend.box = "vertical",
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    axis.text.y = element_text(face = "italic")
  )+
  guides(
    color = guide_legend(order = 1),  # Direction ПЕРВАЯ
    size = guide_legend(order = 2)    # -log10(P-value) ВТОРАЯ
  )

# ==============================================================================
# PART 9: COMBINE & SAVE
# ==============================================================================

combined14 <- (p_abund | p_prev) +
  plot_annotation(
    title = "Bacterial biomarkers associated with endoscopic findings (E3 / E2) in Ulcerative Colitis patients",
    theme = theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5, margin = margin(b = 15))
    )
  )




















# ==============================================================================
# СТРАТИФИЦИРОВАННЫЙ АНАЛИЗ ПО КОМБИНАЦИЯМ ЛЕКАРСТВ
# ==============================================================================

library(maaslin3)
library(dplyr)
library(ggplot2)
library(forcats)
library(tidyr)
library(patchwork)
library(ggrepel)

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
# PART 2: PREPARE METADATA
# ==============================================================================

# Переименовываем лекарства
colnames(input_metadata)[colnames(input_metadata) == "drug_1"] <- "5-ASA"
colnames(input_metadata)[colnames(input_metadata) == "drug_2"] <- "Glucocorticoids"
colnames(input_metadata)[colnames(input_metadata) == "drug_3"] <- "Immunosuppressors"
colnames(input_metadata)[colnames(input_metadata) == "drug_4"] <- "Genetically_Engineered_Therapies"
colnames(input_metadata)[colnames(input_metadata) == "drug_5"] <- "Inhibitor_JAK_kinase"
colnames(input_metadata)[colnames(input_metadata) == "drug_6"] <- "Antibiotics"

# Список лекарств
drugs <- c("5-ASA", "Glucocorticoids", "Immunosuppressors", 
           "Genetically_Engineered_Therapies", "Inhibitor_JAK_kinase", "Antibiotics")

# Преобразуем в numeric
for(drug in drugs) {
  if(drug %in% colnames(input_metadata)) {
    input_metadata[[drug]] <- as.numeric(as.character(input_metadata[[drug]]))
    input_metadata[[drug]][is.na(input_metadata[[drug]])] <- 0
  }
}

# ==============================================================================
# PART 3: ФИЛЬТРУЕМ CONTROL - ОСТАВЛЯЕМ ТОЛЬКО ЗДОРОВЫХ БЕЗ ЛЕКАРСТВ
# ==============================================================================

cat("\n=== BEFORE FILTERING ===\n")
cat("Total samples:", nrow(input_metadata), "\n")
cat("  Control:", sum(input_metadata$diagnosis == "Control"), "\n")
cat("  CD:", sum(input_metadata$diagnosis == "CD"), "\n")
cat("  UC:", sum(input_metadata$diagnosis == "UC"), "\n")

# 1. Удаляем все ND
if("tr" %in% colnames(input_metadata)) {
  nd_samples <- input_metadata$tr == "ND"
  cat("\nND samples:", sum(nd_samples, na.rm = TRUE), "\n")
  cat("  Control ND:", sum(input_metadata$diagnosis == "Control" & nd_samples, na.rm = TRUE), "\n")
  cat("  CD ND:", sum(input_metadata$diagnosis == "CD" & nd_samples, na.rm = TRUE), "\n")
  cat("  UC ND:", sum(input_metadata$diagnosis == "UC" & nd_samples, na.rm = TRUE), "\n")
  
  input_metadata <- input_metadata[!nd_samples, ]
}

# 2. Удаляем Control, которые принимают лекарства
input_metadata <- input_metadata[
  !(input_metadata$diagnosis == "Control" & 
      (input_metadata$`5-ASA` == 1 |
         input_metadata$Glucocorticoids == 1 |
         input_metadata$Immunosuppressors == 1 |
         input_metadata$Genetically_Engineered_Therapies == 1 |
         input_metadata$Inhibitor_JAK_kinase == 1 |
         input_metadata$Antibiotics == 1)),
]

# 3. Удаляем Control с ND (если остались)
if("tr" %in% colnames(input_metadata)) {
  input_metadata <- input_metadata[!(input_metadata$diagnosis == "Control" & input_metadata$tr == "ND"), ]
}

cat("\n=== AFTER FILTERING ===\n")
cat("Total samples:", nrow(input_metadata), "\n")
cat("  Control:", sum(input_metadata$diagnosis == "Control"), "\n")
cat("  CD:", sum(input_metadata$diagnosis == "CD"), "\n")
cat("  UC:", sum(input_metadata$diagnosis == "UC"), "\n")

# Проверяем Control на наличие лекарств
cat("\n=== Control drug usage (should be all 0) ===\n")
control_only <- input_metadata[input_metadata$diagnosis == "Control", ]
for(drug in drugs) {
  if(drug %in% colnames(control_only)) {
    n_drug <- sum(control_only[[drug]] == 1, na.rm = TRUE)
    if(n_drug > 0) {
      cat("  WARNING: Control with", drug, ":", n_drug, "\n")
    } else {
      cat("  ", drug, ": OK (0)\n")
    }
  }
}

# ==============================================================================
# PART 4: ВЫРАВНИВАЕМ ДАННЫЕ
# ==============================================================================

common <- intersect(rownames(input_metadata), rownames(input_data))
input_metadata <- input_metadata[common, , drop = FALSE]
input_data <- input_data[common, , drop = FALSE]

cat("\n=== ALIGNED DATA ===\n")
cat("Samples:", nrow(input_metadata), "\n")
cat("  Control:", sum(input_metadata$diagnosis == "Control"), "\n")
cat("  CD:", sum(input_metadata$diagnosis == "CD"), "\n")
cat("  UC:", sum(input_metadata$diagnosis == "UC"), "\n")

# ==============================================================================
# PART 5: ОПРЕДЕЛЯЕМ ГРУППЫ ДЛЯ СТРАТИФИКАЦИИ
# ==============================================================================

# Создаем переменную для комбинаций
input_metadata$drug_group <- NA

# Для CD
cd_idx <- input_metadata$diagnosis == "CD"
#input_metadata$drug_group[cd_idx & input_metadata$tr == "d100000"] <- "CD_5ASA_mono"
#input_metadata$drug_group[cd_idx & input_metadata$tr == "d000100"] <- "CD_GET_mono"
input_metadata$drug_group[cd_idx & input_metadata$tr == "d001000"] <- "CD_Immunosuppressors_mono"
input_metadata$drug_group[cd_idx & input_metadata$tr == "d011000"] <- "CD_Glucocorticoids_Immunosuppressors"

# Для UC
uc_idx <- input_metadata$diagnosis == "UC"
input_metadata$drug_group[uc_idx & input_metadata$tr == "d100000"] <- "UC_5ASA_mono"
#input_metadata$drug_group[uc_idx & input_metadata$tr == "d110000"] <- "UC_5ASA_Glucocorticoids"
#input_metadata$drug_group[uc_idx & input_metadata$tr == "d111000"] <- "UC_5ASA_Glucocorticoids_Immunosuppressors"

# Control группа (без лекарств)
control_idx <- input_metadata$diagnosis == "Control"
input_metadata$drug_group[control_idx] <- "Control"

# Удаляем NA (пациенты с другими комбинациями)
input_metadata <- input_metadata[!is.na(input_metadata$drug_group), ]

cat("\n=== DRUG GROUPS ===\n")
print(table(input_metadata$drug_group))



# ==============================================================================
# ИСПРАВЛЕННАЯ ФУНКЦИЯ analyze_group
# ==============================================================================

analyze_group <- function(data, metadata, group_name) {
  
  # Берем ВСЕ Control (независимо от drug_group)
  control_data <- metadata[metadata$diagnosis == "Control", ]
  
  # Берем пациентов из нужной группы
  patients_data <- metadata[metadata$drug_group == group_name, ]
  
  # Объединяем
  metadata_sub <- rbind(control_data, patients_data)
  
  # Проверяем, что есть данные
  if(nrow(metadata_sub) == 0) {
    cat("  No data\n")
    return(NULL)
  }
  
  # Создаем переменную для сравнения
  if(grepl("^CD_", group_name)) {
    metadata_sub$comparison <- ifelse(
      metadata_sub$diagnosis == "CD", "CD",
      ifelse(metadata_sub$diagnosis == "Control", "Control", NA)
    )
  } else if(grepl("^UC_", group_name)) {
    metadata_sub$comparison <- ifelse(
      metadata_sub$diagnosis == "UC", "UC",
      ifelse(metadata_sub$diagnosis == "Control", "Control", NA)
    )
  }
  
  metadata_sub <- metadata_sub[!is.na(metadata_sub$comparison), ]
  
  # Проверяем размер
  if(nrow(metadata_sub) < 10) {
    cat("  Not enough samples:", nrow(metadata_sub), "\n")
    return(NULL)
  }
  
  # Проверяем, что есть Control
  n_control <- sum(metadata_sub$comparison == "Control")
  n_patients <- sum(metadata_sub$comparison != "Control")
  
  if(n_control < 3 || n_patients < 3) {
    cat("  Too few Control or patients:", n_control, "/", n_patients, "\n")
    return(NULL)
  }
  
  cat("  Control:", n_control, ", Patients:", n_patients, "\n")
  
  # Выравниваем
  common <- intersect(rownames(metadata_sub), rownames(data))
  metadata_sub <- metadata_sub[common, , drop = FALSE]
  data_sub <- data[common, , drop = FALSE]
  
  # Факторы
  metadata_sub$comparison <- factor(metadata_sub$comparison, levels = c("Control", 
                                                                        ifelse(grepl("^CD_", group_name), "CD", "UC")))
  metadata_sub$sex <- factor(metadata_sub$sex, levels = c("F", "M"))
  metadata_sub$patient <- as.character(metadata_sub$patient)
  metadata_sub$age <- as.numeric(as.character(metadata_sub$age))
  
  # Запускаем MaAsLin3
  output_path <- paste0("Desktop/WORK/gut/1_stage/R/stratified_groups/", group_name)
  
  fit <- maaslin3(
    input_data = data_sub,
    input_metadata = metadata_sub,
    output = output_path,
    
    fixed_effects = c("comparison", "age", "sex"),
    random_effects = c("patient"),
    
    reference = c("comparison=Control", "sex=F"),
    
    normalization = "TSS",
    transform = "LOG",
    correction = "BH",
    
    min_abundance = 0.0001,
    min_prevalence = 0.1,
    small_random_effects = TRUE
  )
  
  # Загружаем результаты
  results_file <- paste0(output_path, "/all_results.tsv")
  
  if(file.exists(results_file)) {
    results <- read.delim(results_file, check.names = FALSE)
    sig <- results %>%
      filter(pval_joint < 0.05) %>%
      filter(metadata == "comparison")
    
    return(sig)
  }
  
  return(NULL)
}

# ==============================================================================
# ПРОВЕРКА: Сколько Control доступно
# ==============================================================================

cat("\n=== CONTROL SAMPLES ===\n")
control_data <- input_metadata[input_metadata$diagnosis == "Control", ]
cat("Total Control:", nrow(control_data), "\n")

# Проверяем, что Control не имеют лекарств
for(drug in drugs) {
  if(drug %in% colnames(control_data)) {
    n_drug <- sum(control_data[[drug]] == 1, na.rm = TRUE)
    if(n_drug > 0) {
      cat("  WARNING: Control with", drug, ":", n_drug, "\n")
    }
  }
}

# ==============================================================================
# ЗАПУСК АНАЛИЗА ДЛЯ ВСЕХ ГРУПП
# ==============================================================================

dir.create("Desktop/WORK/gut/1_stage/R/stratified_groups", showWarnings = FALSE)

groups <- c("CD_Immunosuppressors_mono", 
            "CD_Glucocorticoids_Immunosuppressors",
            "UC_5ASA_mono")

all_results <- list()

cat("\n========================================\n")
cat("RUNNING STRATIFIED ANALYSIS\n")
cat("========================================\n")

for(group in groups) {
  cat("\n", group, ":\n")
  res <- analyze_group(input_data, input_metadata, group)
  all_results[[group]] <- res
  
  if(!is.null(res)) {
    cat("  Significant associations:", nrow(res), "\n")
    if(nrow(res) > 0) {
      print(res[, c("feature", "model", "coef", "pval_joint")])
    }
  } else {
    cat("  No results\n")
  }
}





# ==============================================================================
# ВСЕ ГРУППЫ В ОДНОМ ГРАФИКЕ
# ==============================================================================

# Объединяем все значимые результаты из всех групп
all_significant <- data.frame()

for(group in names(all_results)) {
  if(!is.null(all_results[[group]]) && nrow(all_results[[group]]) > 0) {
    res <- all_results[[group]]
    res$group <- group
    all_significant <- rbind(all_significant, res)
  }
}

if(nrow(all_significant) > 0) {
  
  # Красивые названия групп
  group_labels <- c(
    "CD_Immunosuppressors_mono" = "CD + Immunosuppressors",
    "CD_Glucocorticoids_Immunosuppressors" = "CD + Gluco + Immuno",
    "UC_5ASA_mono" = "UC + 5-ASA"
  )
  
  all_significant$group_label <- group_labels[all_significant$group]
  all_significant$neglogp <- -log10(all_significant$pval_joint)
  all_significant$direction <- ifelse(all_significant$coef > 0, "Higher in patients", "Higher in Control")
  
  # График с фасетами
  p_all <- ggplot(all_significant, aes(x = coef, y = feature, color = direction, size = neglogp)) +
    geom_point(alpha = 0.8) +
    geom_errorbarh(aes(xmin = coef - 1.96*0.3, xmax = coef + 1.96*0.3), 
                   height = 0.2, alpha = 0.5) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    facet_grid(model ~ group_label, scales = "free_y") +
    scale_color_manual(
      name = "Direction",
      values = c("Higher in patients" = "#1B4F72", "Higher in Control" = "#C06014")
    ) +
    scale_size_continuous(name = "-log10(P-value)", range = c(2, 6)) +
    scale_y_discrete(labels = function(x) gsub("_", " ", x)) +
    labs(
      title = "Bacterial biomarkers stratified by therapy",
      subtitle = "Significant associations (p < 0.05)",
      x = "β coefficient (±95% CI)",
      y = ""
    ) +
    theme_minimal() +
    theme(
      strip.text = element_text(face = "bold", size = 11),
      axis.text.y = element_text(face = "italic", size = 9),
      legend.position = "bottom",
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      plot.subtitle = element_text(hjust = 0.5, size = 11),
      panel.grid.minor = element_blank()
    )
  
  print(p_all)
  
  ggsave(
    "Desktop/WORK/gut/1_stage/R/stratified_groups/all_groups_results.pdf",
    p_all,
    width = 14, height = 10, device = cairo_pdf
  )
}













###########ВСЕ ГРАФИКИ РАЗДЕЛЬНО

cd_immuno_results <- all_significant %>%
  filter(group == "CD_Immunosuppressors_mono") %>% filter(pval_joint < 0.05)

# ==============================================================================
# PART 7: PLOT (ABUNDANCE) - Activity
# ==============================================================================

abund <- cd_immuno_results %>%
  filter(model == "abundance") %>%
  filter(!is.na(coef), !is.na(stderr), !is.na(pval_joint)) %>%
  mutate(
    feature = fct_reorder(feature, coef),
    direction = ifelse(coef > 0, "Higher in Immunosuppressors", "Higher in Control"),
    lower = coef - 1.96 * stderr,
    upper = coef + 1.96 * stderr,
    neglogp = -log10(pval_joint)
  )

p_abund <- ggplot(abund, aes(x = coef, y = feature, color = direction)) +
  geom_point(aes(size = neglogp), alpha = 0.8) +
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.2, alpha = 0.6) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_color_manual(
    name = "Direction",
    values = c("Higher in Immunosuppressors" = "#1B4F72", "Higher in Control" = "#C06014")
  ) +
  scale_size_continuous(
    name = "-log10(P-value)"
  ) +
  scale_y_discrete(
    labels = function(x) gsub("_", " ", x)
  ) +
  labs(
    title = "Abundance associations",
    subtitle = "CD patients (Immunosuppressors) vs Control\nadjusted for age and gender",
    x = "β coefficient (±95% CI)", 
    y = ""
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    legend.box = "vertical",
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    axis.text.y = element_text(face = "italic")
  )+
  guides(
    color = guide_legend(order = 1),  # Direction ПЕРВАЯ
    size = guide_legend(order = 2)    # -log10(P-value) ВТОРАЯ
  )

# ==============================================================================
# PART 8: PLOT (PREVALENCE)
# ==============================================================================

prev <- cd_immuno_results %>%
  filter(model == "prevalence") %>%
  mutate(
    feature = fct_reorder(feature, coef),
    direction = ifelse(coef > 0, "More frequent in Immunosuppressors", "More frequent in Control"),
    lower = coef - 1.96 * stderr,
    upper = coef + 1.96 * stderr,
    neglogp = -log10(pval_joint)
  )

p_prev <- ggplot(prev, aes(x = coef, y = feature, color = direction)) +
  geom_point(aes(size = neglogp), alpha = 0.8) +
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.2, alpha = 0.6) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_color_manual(
    name = "Direction",
    values = c("More frequent in Immunosuppressors" = "#2E8B57", "More frequent in Control" = "#B22222")
  ) +
  scale_size_continuous(
    name = "-log10(P-value)"
  ) +
  scale_y_discrete(
    labels = function(x) gsub("_", " ", x)
  ) +
  labs(
    title = "Prevalence associations",
    subtitle = "CD patients (Immunosuppressors) vs Control\nadjusted for age and gender",
    x = "β coefficient (±95% CI)", 
    y = ""
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    legend.box = "vertical",
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    axis.text.y = element_text(face = "italic")
  )+
  guides(
    color = guide_legend(order = 1),  # Direction ПЕРВАЯ
    size = guide_legend(order = 2)    # -log10(P-value) ВТОРАЯ
  )

# ==============================================================================
# PART 9: COMBINE & SAVE
# ==============================================================================

combined15 <- (p_abund | p_prev) +
  plot_annotation(
    title = "Bacterial biomarkers associated with Crohn's disease patients (Immunosuppressors) vs Control",
    theme = theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5, margin = margin(b = 15))
    )
  )




###########ВСЕ ГРАФИКИ РАЗДЕЛЬНО

cd_glu_immuno_results <- all_significant %>%
  filter(group == "CD_Glucocorticoids_Immunosuppressors") %>% filter(pval_joint < 0.05)

# ==============================================================================
# PART 7: PLOT (ABUNDANCE) - Activity
# ==============================================================================

abund <- cd_glu_immuno_results %>%
  filter(model == "abundance") %>%
  filter(!is.na(coef), !is.na(stderr), !is.na(pval_joint)) %>%
  mutate(
    feature = fct_reorder(feature, coef),
    direction = ifelse(coef > 0, "Higher in G + I", "Higher in Control"),
    lower = coef - 1.96 * stderr,
    upper = coef + 1.96 * stderr,
    neglogp = -log10(pval_joint)
  )

p_abund <- ggplot(abund, aes(x = coef, y = feature, color = direction)) +
  geom_point(aes(size = neglogp), alpha = 0.8) +
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.2, alpha = 0.6) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_color_manual(
    name = "Direction",
    values = c("Higher in G + I" = "#1B4F72", "Higher in Control" = "#C06014")
  ) +
  scale_size_continuous(
    name = "-log10(P-value)"
  ) +
  scale_y_discrete(
    labels = function(x) gsub("_", " ", x)
  ) +
  labs(
    title = "Abundance associations",
    subtitle = "CD patients\n(Glucocorticoids + Immunosuppressors) vs Control\nadjusted for age and gender",
    x = "β coefficient (±95% CI)", 
    y = ""
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    legend.box = "vertical",
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    axis.text.y = element_text(face = "italic")
  )+
  guides(
    color = guide_legend(order = 1),  # Direction ПЕРВАЯ
    size = guide_legend(order = 2)    # -log10(P-value) ВТОРАЯ
  )

# ==============================================================================
# PART 8: PLOT (PREVALENCE)
# ==============================================================================

prev <- cd_glu_immuno_results %>%
  filter(model == "prevalence") %>%
  mutate(
    feature = fct_reorder(feature, coef),
    direction = ifelse(coef > 0, "More frequent in G + I", "More frequent in Control"),
    lower = coef - 1.96 * stderr,
    upper = coef + 1.96 * stderr,
    neglogp = -log10(pval_joint)
  )

p_prev <- ggplot(prev, aes(x = coef, y = feature, color = direction)) +
  geom_point(aes(size = neglogp), alpha = 0.8) +
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.2, alpha = 0.6) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_color_manual(
    name = "Direction",
    values = c("More frequent in G + I" = "#2E8B57", "More frequent in Control" = "#B22222")
  ) +
  scale_size_continuous(
    name = "-log10(P-value)"
  ) +
  scale_y_discrete(
    labels = function(x) gsub("_", " ", x)
  ) +
  labs(
    title = "Prevalence associations",
    subtitle = "CD patients\n(Glucocorticoids + Immunosuppressors) vs Control\nadjusted for age and gender",
    x = "β coefficient (±95% CI)", 
    y = ""
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    legend.box = "vertical",
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    axis.text.y = element_text(face = "italic")
  )+
  guides(
    color = guide_legend(order = 1),  # Direction ПЕРВАЯ
    size = guide_legend(order = 2)    # -log10(P-value) ВТОРАЯ
  )

# ==============================================================================
# PART 9: COMBINE & SAVE
# ==============================================================================

combined16 <- (p_abund | p_prev) +
  plot_annotation(
    title = "Bacterial biomarkers associated with\nCrohn's disease patients (Glucocorticoids + Immunosuppressors) vs Control",
    theme = theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5, margin = margin(b = 15))
    )
  )





###########ВСЕ ГРАФИКИ РАЗДЕЛЬНО

uc_5asa_results <- all_significant %>%
  filter(group == "UC_5ASA_mono") %>% filter(pval_joint < 0.05)

# ==============================================================================
# PART 7: PLOT (ABUNDANCE) - Activity
# ==============================================================================

abund <- uc_5asa_results %>%
  filter(model == "abundance") %>%
  filter(!is.na(coef), !is.na(stderr), !is.na(pval_joint)) %>%
  mutate(
    feature = fct_reorder(feature, coef),
    direction = ifelse(coef > 0, "Higher in 5-ASA", "Higher in Control"),
    lower = coef - 1.96 * stderr,
    upper = coef + 1.96 * stderr,
    neglogp = -log10(pval_joint)
  )

p_abund <- ggplot(abund, aes(x = coef, y = feature, color = direction)) +
  geom_point(aes(size = neglogp), alpha = 0.8) +
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.2, alpha = 0.6) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_color_manual(
    name = "Direction",
    values = c("Higher in 5-ASA" = "#1B4F72", "Higher in Control" = "#C06014")
  ) +
  scale_size_continuous(
    name = "-log10(P-value)"
  ) +
  scale_y_discrete(
    labels = function(x) gsub("_", " ", x)
  ) +
  labs(
    title = "Abundance associations",
    subtitle = "UC patients (5-ASA) vs Control\nadjusted for age and gender",
    x = "β coefficient (±95% CI)", 
    y = ""
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    legend.box = "vertical",
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    axis.text.y = element_text(face = "italic")
  )+
  guides(
    color = guide_legend(order = 1),  # Direction ПЕРВАЯ
    size = guide_legend(order = 2)    # -log10(P-value) ВТОРАЯ
  )

# ==============================================================================
# PART 8: PLOT (PREVALENCE)
# ==============================================================================

prev <- uc_5asa_results %>%
  filter(model == "prevalence") %>%
  mutate(
    feature = fct_reorder(feature, coef),
    direction = ifelse(coef > 0, "More frequent in 5-ASA", "More frequent in Control"),
    lower = coef - 1.96 * stderr,
    upper = coef + 1.96 * stderr,
    neglogp = -log10(pval_joint)
  )

p_prev <- ggplot(prev, aes(x = coef, y = feature, color = direction)) +
  geom_point(aes(size = neglogp), alpha = 0.8) +
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.2, alpha = 0.6) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_color_manual(
    name = "Direction",
    values = c("More frequent in 5-ASA" = "#2E8B57", "More frequent in Control" = "#B22222")
  ) +
  scale_size_continuous(
    name = "-log10(P-value)"
  ) +
  scale_y_discrete(
    labels = function(x) gsub("_", " ", x)
  ) +
  labs(
    title = "Prevalence associations",
    subtitle = "UC patients (5-ASA) vs Control\nadjusted for age and gender",
    x = "β coefficient (±95% CI)", 
    y = ""
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    legend.box = "vertical",
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    axis.text.y = element_text(face = "italic")
  )+
  guides(
    color = guide_legend(order = 1), 
    size = guide_legend(order = 2) 
  )

# ==============================================================================
# PART 9: COMBINE & SAVE
# ==============================================================================

combined17 <- (p_abund | p_prev) +
  plot_annotation(
    title = "Bacterial biomarkers associated with Ulcerative Colitis patients (5-ASA) vs Control",
    theme = theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5, margin = margin(b = 15))
    )
  )







combined_all <- wrap_elements(combined1) / wrap_elements(combined2) / wrap_elements(combined3) / wrap_elements(combined4) / wrap_elements(combined5) / wrap_elements(combined7) / wrap_elements(combined8) / wrap_elements(combined9) / wrap_elements(combined10) / wrap_elements(combined11) / wrap_elements(combined13) / wrap_elements(combined14) / wrap_elements(combined15) / wrap_elements(combined16) / wrap_elements(combined17)

ggsave(
  "Desktop/WORK/gut/1_stage/R/biomarkers.pdf",
  combined_all,
  width = 12, height = 75, device = cairo_pdf,
  limitsize = FALSE
)
