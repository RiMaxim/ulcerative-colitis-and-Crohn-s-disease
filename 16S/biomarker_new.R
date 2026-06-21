####################PART1


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
colnames(metadata_filtered)[colnames(metadata_filtered) == "drug_1"] <- "5-ASA"
colnames(metadata_filtered)[colnames(metadata_filtered) == "drug_2"] <- "Glucocorticoids"
colnames(metadata_filtered)[colnames(metadata_filtered) == "drug_3"] <- "Immunosuppressors"
colnames(metadata_filtered)[colnames(metadata_filtered) == "drug_4"] <- "Genetically_Engineered_Therapies"
colnames(metadata_filtered)[colnames(metadata_filtered) == "drug_5"] <- "Inhibitor_JAK_kinase"
colnames(metadata_filtered)[colnames(metadata_filtered) == "drug_6"] <- "Antibiotics"

# ensure numeric (important for MaAsLin3)
drug_cols <- c("5-ASA","Glucocorticoids","Immunosuppressors","Genetically_Engineered_Therapies","Inhibitor_JAK_kinase","Antibiotics")
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
    direction = ifelse(coef > 0, "Higher in CD", "Lower in CD"),
    lower = coef - 1.96 * stderr,
    upper = coef + 1.96 * stderr,
    neglogp = -log10(pval_joint)
  )

p_abund <- ggplot(abund, aes(x = coef, y = feature, color = direction)) +
  geom_point(aes(size = neglogp), alpha = 0.8) +
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.2, alpha = 0.6) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_color_manual(
  name = "Direction",  # ЯВНО УКАЗЫВАЕМ НАЗВАНИЕ
  values = c("Higher in CD" = "#1B4F72", "Lower in CD" = "#C06014")
) +
  scale_size_continuous(
    name = "-log10(P-value)",  # ЯВНО УКАЗЫВАЕМ НАЗВАНИЕ
  ) +
  scale_y_discrete(
    labels = function(x) gsub("_", " ", x)  # ЗАМЕНЯЕМ _ НА ПРОБЕЛ
  ) +
  labs(
    title = "Abundance associations",
    subtitle = "CD vs Control\nadjusted for age, sex",
    x = "β coefficient (±95% CI)", y = ""
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    legend.box = "vertical",
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    axis.text.y = element_text(face = "italic")  # КУРСИВ
  )



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
  scale_color_manual(
  name = "Direction",  # ТО ЖЕ НАЗВАНИЕ
  values = c("More frequent in CD" = "#2E8B57", "Less frequent in CD" = "#B22222")
) +
  scale_size_continuous(
    name = "-log10(P-value)",  # ТО ЖЕ НАЗВАНИЕ
  ) +
  scale_y_discrete(
    labels = function(x) gsub("_", " ", x)  # ЗАМЕНЯЕМ _ НА ПРОБЕЛ
  ) +
  labs(
    title = "Prevalence associations",
    subtitle = "CD vs Control\nadjusted for age, sex",
    x = "β coefficient (±95% CI)", y = ""
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    legend.box = "vertical",
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    axis.text.y = element_text(face = "italic")  # КУРСИВ
  )


# ==============================================================================
# PART 9: COMBINE & SAVE
# ==============================================================================

combined1 <- (p_abund | p_prev) +
  plot_annotation(
    title = "Bacterial biomarkers associated with Crohn's disease",
    theme = theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5, margin = margin(b = 15)))
  )




################PART2



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
    direction = ifelse(coef > 0, "Higher in UC", "Lower in UC"),
    lower = coef - 1.96 * stderr,
    upper = coef + 1.96 * stderr,
    neglogp = -log10(pval_joint)
  )

p_abund <- ggplot(abund, aes(x = coef, y = feature, color = direction)) +
  geom_point(aes(size = neglogp), alpha = 0.8) +
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.2, alpha = 0.6) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_color_manual(
    name = "Direction",  # ЯВНО УКАЗЫВАЕМ НАЗВАНИЕ
    values = c("Higher in UC" = "#1B4F72", "Lower in UC" = "#C06014")
  ) +
  scale_size_continuous(
    name = "-log10(P-value)",  # ЯВНО УКАЗЫВАЕМ НАЗВАНИЕ
  ) +
  scale_y_discrete(
    labels = function(x) gsub("_", " ", x)  # ЗАМЕНЯЕМ _ НА ПРОБЕЛ
  ) +
  labs(
    title = "Abundance associations",
    subtitle = "UC vs Control\nadjusted for age, sex",
    x = "β coefficient (±95% CI)", y = ""
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    legend.box = "vertical",
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    axis.text.y = element_text(face = "italic")  # КУРСИВ
  )

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
  scale_color_manual(
    name = "Direction",  # ТО ЖЕ НАЗВАНИЕ
    values = c("More frequent in UC" = "#2E8B57", "Less frequent in UC" = "#B22222")
  ) +
  scale_size_continuous(
    name = "-log10(P-value)",  # ТО ЖЕ НАЗВАНИЕ
  ) +
  scale_y_discrete(
    labels = function(x) gsub("_", " ", x)  # ЗАМЕНЯЕМ _ НА ПРОБЕЛ
  ) +
  labs(
    title = "Prevalence associations",
    subtitle = "UC vs Control\nadjusted for age, sex",
    x = "β coefficient (±95% CI)", y = ""
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    legend.box = "vertical",
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    axis.text.y = element_text(face = "italic")  # КУРСИВ
  )

# ==============================================================================
# PART 9: SAVE
# ==============================================================================

combined2 <- (p_abund | p_prev) +
  plot_annotation(
    title = "Bacterial biomarkers associated with Ulcerative Colitis",
    theme = theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5, margin = margin(b = 15)))
  )







######################PART3



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
dim(metadata_filtered)
metadata_filtered <- metadata_filtered[complete.cases(metadata_filtered[drug_cols]), ]
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
    subtitle = "CD vs UC\nadjusted for age, sex, therapy and extraintestinal manifistation",
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
    subtitle = "CD vs UC\nadjusted for age, sex, therapy and extraintestinal manifistation",
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
dim(metadata_filtered)


# Преобразуем er и em в факторы
metadata_filtered$er <- factor(metadata_filtered$er, levels = c(0, 1))
metadata_filtered$em <- factor(metadata_filtered$em, levels = c(0, 1))

# Создаем комбинированную переменную er_em для взаимодействия
# или используем er с поправкой на em

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
    subtitle = "CD patients only: exacerbation vs remission\nadjusted for age, sex, therapy and extraintestinal manifistation",
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
    subtitle = "CD patients only: exacerbation vs remission\nadjusted for age, sex, therapy and extraintestinal manifistation",
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
# PART 2: FILTER UC PATIENTS AND CREATE ER VARIABLE
# ==============================================================================

# Оставляем только UC
metadata_filtered <- input_metadata[input_metadata$diagnosis == "UC", ]

# Удаляем образцы с ER=1 & EM=0 (только 1 образец)
samples_to_remove <- rownames(metadata_filtered)[metadata_filtered$er == 1 & metadata_filtered$em == 0]
if(length(samples_to_remove) > 0) {
  cat("Removing", length(samples_to_remove), "samples with ER=1 & EM=0\n")
  metadata_filtered <- metadata_filtered[!rownames(metadata_filtered) %in% samples_to_remove, ]
  data_filtered <- data_filtered[!rownames(data_filtered) %in% samples_to_remove, ]
}

# Преобразуем er и em в факторы
metadata_filtered$er <- factor(metadata_filtered$er, levels = c(0, 1))
metadata_filtered$em <- factor(metadata_filtered$em, levels = c(0, 1))

# Создаем комбинированную переменную er_em для взаимодействия
# или используем er с поправкой на em

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

# ==============================================================================
# PART 4: ALIGN DATA
# ==============================================================================

common_samples <- intersect(rownames(metadata_filtered), rownames(input_data))

metadata_filtered <- metadata_filtered[common_samples, , drop = FALSE]
data_filtered <- input_data[common_samples, , drop = FALSE]

metadata_filtered$sex <- factor(metadata_filtered$sex, levels = c("F", "M"))

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
    subtitle = "UC patients only: exacerbation vs remission\nadjusted for age, sex, therapy and extraintestinal manifistation",
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
    subtitle = "UC patients only: exacerbation vs remission\nadjusted for age, sex, therapy and extraintestinal manifistation",
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

combined5 <- (p_abund | p_prev) +
  plot_annotation(
    title = "Bacterial biomarkers associated with exacerbation / remission status in Ulcerative Colitis patients",
    theme = theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5, margin = margin(b = 15))
    )
  )



combined_all <- wrap_elements(combined1) / wrap_elements(combined2) / wrap_elements(combined3) / wrap_elements(combined4) / wrap_elements(combined5)

ggsave(
  "Desktop/WORK/gut/1_stage/R/biomarkers.pdf",
  combined_all,
  width = 12, height = 25, device = cairo_pdf
)
