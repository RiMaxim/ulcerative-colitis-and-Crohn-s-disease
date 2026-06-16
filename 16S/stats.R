# ============================================================================
# Script for visualizing read quality before and after filtering
# Input file: table with statistics per sample
# ============================================================================

# Load required libraries
library(ggplot2)
library(gridExtra)
library(dplyr)

# ============================================================================
# 1. Read input file
# ============================================================================

# Specify the path to your input file
input_file <- "Desktop/WORK/gut/1_stage/R/stats.tsv"  # <-- CHANGE THIS PATH


# Expected structure of input file (TSV with headers):
# SampleID | Reads_raw | Reads_filtered | Quality_raw | Quality_filtered | Length_raw | Length_filtered | SD_raw | SD_filtered
#
# Example:
# Sample1 | 45000 | 42000 | 28.5 | 30.2 | 1200 | 1150 | 450 | 420
# Sample2 | 12000 | 11000 | 27.8 | 29.5 | 1180 | 1130 | 430 | 400

data <- read.table(input_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Check data structure
cat("Input file structure:\n")
str(data)

# Check for required columns
required_cols <- c("SampleID", "Reads_raw", "Reads_filtered", "Quality_raw", 
                   "Quality_filtered", "Length_raw", "Length_filtered", 
                   "SD_raw", "SD_filtered")

missing_cols <- setdiff(required_cols, colnames(data))
if (length(missing_cols) > 0) {
  stop("Missing columns in file: ", paste(missing_cols, collapse = ", "))
}

# ============================================================================
# 2. Generate plots
# ============================================================================

# 2.1 Number of reads before and after filtering
p1 <- ggplot(data) +
  geom_segment(aes(x = reorder(SampleID, Reads_raw), 
                   xend = reorder(SampleID, Reads_raw),
                   y = Reads_raw, yend = Reads_filtered),
               color = "gray50", size = 0.5) +
  geom_hline(yintercept = 10000, color = "red", linetype = "dashed", size = 0.8) +
  geom_hline(yintercept = 5000, color = "pink", linetype = "dashed", size = 0.8) +
  geom_point(aes(x = reorder(SampleID, Reads_raw), y = Reads_raw, 
                 color = "Before filtering (raw)"), 
             size = 2) +
  geom_point(aes(x = reorder(SampleID, Reads_raw), y = Reads_filtered, 
                 color = "After filtering (filtered)"), 
             size = 2) +
  scale_color_manual(values = c(
    "Before filtering (raw)" = "darkred",
    "After filtering (filtered)" = "steelblue"
  )) +
  labs(title = "Number of reads before and after filtering",
       x = "Sample", y = "Number of reads",
       color = "Stage") +
  scale_y_continuous(labels = scales::comma) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 3),
        legend.position = "bottom",
        plot.title = element_text(face = "bold", hjust = 0.5),
        axis.title   = element_text(size = 16),
        axis.text    = element_text(size = 12))+ylim(0,400000)

# 2.2 Read quality plot (Q-score)
p2 <- ggplot(data) +
  geom_segment(aes(x = reorder(SampleID, Quality_raw), 
                   xend = reorder(SampleID, Quality_raw),
                   y = Quality_raw, yend = Quality_filtered),
               color = "gray70", size = 0.5) +
  geom_point(aes(x = reorder(SampleID, Quality_raw), y = Quality_raw, 
                 color = "Before filtering"), size = 2) +
  geom_point(aes(x = reorder(SampleID, Quality_raw), y = Quality_filtered, 
                 color = "After filtering"), size = 2) +
  scale_color_manual(values = c("Before filtering" = "darkred", 
                                "After filtering" = "steelblue")) +
  labs(title = "Read quality (Q-score)",
       x = "Sample", y = "Mean quality", 
       color = "Stage") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 3),
        legend.position = "bottom",
        plot.title = element_text(face = "bold", hjust = 0.5),
        axis.title   = element_text(size = 16),
        axis.text    = element_text(size = 12))

# 2.3 Read length plot
p3 <- ggplot(data) +
  geom_segment(aes(x = reorder(SampleID, Length_raw), 
                   xend = reorder(SampleID, Length_raw),
                   y = Length_raw, yend = Length_filtered),
               color = "gray70", size = 0.5) +
  geom_point(aes(x = reorder(SampleID, Length_raw), y = Length_raw, 
                 color = "Before filtering"), size = 2) +
  geom_point(aes(x = reorder(SampleID, Length_raw), y = Length_filtered, 
                 color = "After filtering"), size = 2) +
  scale_color_manual(values = c("Before filtering" = "darkred", 
                                "After filtering" = "steelblue")) +
  labs(title = "Read length",
       x = "Sample", y = "Mean length (bp)",
       color = "Stage") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 3),
        legend.position = "bottom",
        plot.title = element_text(face = "bold", hjust = 0.5),
        axis.title   = element_text(size = 16),
        axis.text    = element_text(size = 12))

# 2.4 Read length variability plot (standard deviation)
p4 <- ggplot(data) +
  geom_segment(aes(x = reorder(SampleID, SD_raw), 
                   xend = reorder(SampleID, SD_raw),
                   y = SD_raw, yend = SD_filtered),
               color = "gray70", size = 0.5) +
  geom_point(aes(x = reorder(SampleID, SD_raw), y = SD_raw, 
                 color = "Before filtering"), size = 2) +
  geom_point(aes(x = reorder(SampleID, SD_raw), y = SD_filtered, 
                 color = "After filtering"), size = 2) +
  scale_color_manual(values = c("Before filtering" = "darkred", 
                                "After filtering" = "steelblue")) +
  labs(title = "Read length variability",
       x = "Sample", y = "Standard deviation",
       color = "Stage") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 3),
        legend.position = "bottom",
        plot.title = element_text(face = "bold", hjust = 0.5),
        axis.title   = element_text(size = 16),
        axis.text    = element_text(size = 12))

# ============================================================================
# 3. Combine and save plots
# ============================================================================

# Combine all plots into a 2x2 grid
combined_plot <- grid.arrange(p1, p2, p3, p4, ncol = 2, nrow = 2)

# Create output directory if it doesn't exist
output_dir <- "Desktop/WORK/gut/1_stage/R/"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Save combined plot as PDF
output_file <- file.path(output_dir, "192_files_reads_quality_length_SD.pdf")
ggsave(output_file, combined_plot, width = 16, height = 10)

cat("Plot saved to:", output_file, "\n")

