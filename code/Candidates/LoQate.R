# =============================================================================
# Author:            Jaime Martínez Cazón
#
# Description:
# This script analyzes changes in protein expression noise in S. cerevisiae
# across different environmental stress conditions (DTT, H2O2, and nutrient
# starvation), using data derived from the LoQate resource.
#
# The script performs the following key tasks:
# 1. Data Loading and Pre-processing: It loads datasets containing protein
#    expression statistics (median, STD, fold change) for control and
#    stress conditions. It includes a commented-out section detailing the
#    initial data filtering for cytoplasmic proteins and the calculation of
#    normalized noise metrics.
# 2. Outlier Identification: It uses a custom function to identify outlier
#    proteins in scatter plots. Outliers are defined as proteins with the
#    largest Euclidean distance from the center of the data cloud (top 1%).
# 3. Visualization: It generates two types of scatter plots for each stress
#    condition to investigate how stress impacts protein expression:
#    a) Baseline Noise vs. Fold Change: To see if a protein's inherent
#       noise level predicts the magnitude of its expression change under stress.
#    b) Baseline Noise vs. Stress Noise: To directly compare noise levels
#       before and after stress application.
#
# =============================================================================


# =============================================================================
# LOAD LIBRARIES
# =============================================================================
library(ggplot2)      # For creating advanced plots and visualizations.
library(dplyr)        # For data manipulation and transformation.
library(stringr)      # For string and pattern matching operations.


# =============================================================================
# SECTION 1: DATA LOADING AND PREPARATION
# =============================================================================

# --- 1.1 Load Pre-processed Data Files ---
# These .rds files contain pre-filtered and processed noise data for each condition.
# The generation of these files is detailed in the commented-out block below.
# Paths are relative to the project's root directory.
noise_DTT <- readRDS("data/Different conditions/noise_DTT.rds")
noise_H2O2 <- readRDS("data/Different conditions/noise_H2O2.rds")
noise_STARVATION <- readRDS("data/Different conditions/noise_STARVATION.rds")


# --- 1.2 (Commented Out) Initial Data Processing Workflow ---
# This block shows how the .rds files were created from the raw LoQate data.
# It is kept for reproducibility but is not executed in a standard run.
'
# --- Load raw annotation and conditions data ---
protein_info <- read.csv("data/Single cell/summarySC.csv")
Conditions_orf <- read.csv("data/Different conditions/NoiseData_DifferentConditions.csv")

# --- Filter for cytoplasmic proteins ---
# This step is crucial for data quality, retaining proteins localized to the
# cytoplasm where 2D fluorescence measurements are most reliable.
cytoplasm_list <- protein_info[
  grepl("cytoplasm", protein_info$SUBCELLULAR.LOCATIONS, ignore.case = TRUE) |
    grepl("cytosol", protein_info$SUBCELLULAR.LOCATIONS, ignore.case = TRUE) |
    is.na(protein_info$SUBCELLULAR.LOCATIONS), 
]
cytoplasm_orfs <- unique(cytoplasm_list$name)

# Apply the cytoplasmic filter to the main conditions dataset.
Conditions_orf_cyto <- Conditions_orf %>%
  filter(ORF %in% cytoplasm_orfs)

# --- Calculate Normalized Noise ---
# Calculate the coefficient of variation (STD / Median) for each condition.
# This normalized metric allows for direct comparison of noise across proteins
# with different absolute expression levels.
Conditions_orf_cyto <- Conditions_orf_cyto %>%
  mutate(
    Control.STD.Normalized = Control.STD / Control.Median,
    DTT.STD.Normalized = DTT.STD / DTT.Median,
    H2O2.STD.Normalized = H2O2.STD / H2O2.Median,
    Starvation.STD.Normalized = Starvation.STD / Starvation.Median
  )

# --- Create and Save Condition-Specific DataFrames ---
noise_DTT <- Conditions_orf_cyto %>%
  select(ORF, Gene, Control.STD.Normalized, DTT.STD.Normalized, DTT.FoldChange)

noise_H2O2 <- Conditions_orf_cyto %>%
  select(ORF, Gene, Control.STD.Normalized, H2O2.STD.Normalized, H2O2.FoldChange)

noise_STARVATION <- Conditions_orf_cyto %>%
  select(ORF, Gene, Control.STD.Normalized, Starvation.STD.Normalized, Starvation.FoldChange)

# Save the final data objects for easy loading in future sessions.
saveRDS(noise_DTT, "data/Different conditions/noise_DTT.rds")
saveRDS(noise_H2O2, "data/Different conditions/noise_H2O2.rds")
saveRDS(noise_STARVATION, "data/Different conditions/noise_STARVATION.rds")
'


# =============================================================================
# SECTION 2: UTILITY FUNCTIONS FOR OUTLIER DETECTION
# =============================================================================

#' @description Calculates the Euclidean distance of a point (x, y) from a center point.
calculate_distance <- function(x, y, center_x, center_y) {
  sqrt((x - center_x)^2 + (y - center_y)^2)
}

#' @description Identifies outlier data points in a 2D plot and prepares them for highlighting.
#' @param df The input data frame.
#' @param y_col The name of the column to be used for the y-axis. This allows the
#'        function to be used flexibly for both noise vs. noise and noise vs. fold change plots.
#' @return The input data frame with two new columns: 'distance' and 'color',
#'         where outliers (top 1% by distance) are marked 'red'.
highlight_outliers <- function(df, y_col) {
  # Calculate the center of the data cloud.
  center_x <- mean(df$Control.STD.Normalized, na.rm = TRUE)
  center_y <- mean(df[[y_col]], na.rm = TRUE)
  
  # Calculate the distance of each point from the center.
  df$distance <- mapply(calculate_distance, df$Control.STD.Normalized, df[[y_col]],
                        MoreArgs = list(center_x = center_x, center_y = center_y))
  
  # Define the outlier threshold as the 99th percentile of distances.
  distance_threshold <- quantile(df$distance, 0.99, na.rm = TRUE)
  
  # Assign color based on whether the point's distance exceeds the threshold.
  df$color <- ifelse(df$distance > distance_threshold, "red", "black")
  
  return(df)
}


# =============================================================================
# SECTION 3: VISUALIZATION I - NOISE vs. FOLD CHANGE
# =============================================================================
# These plots investigate whether a protein's baseline noise level in control
# conditions is related to the magnitude of its expression change under stress.
# Each point is a protein, labeled by its gene name.

# --- DTT Stress ---
ggplot(highlight_outliers(noise_DTT, "DTT.FoldChange"), aes(x = DTT.FoldChange, y = Control.STD.Normalized, label = Gene)) +
  geom_text(aes(color = color), fontface = "bold", size = 3.5, show.legend = FALSE) +
  labs(title = "DTT Stress", x = "Median Fold Change", y = "Control Noise (Normalized STD)") +
  scale_color_manual(values = c("black" = "black", "red" = "#E41A1C")) +
  theme_minimal(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16))

# --- H2O2 Stress ---
ggplot(highlight_outliers(noise_H2O2, "H2O2.FoldChange"), aes(x = H2O2.FoldChange, y = Control.STD.Normalized, label = Gene)) +
  geom_text(aes(color = color), fontface = "bold", size = 3.5, show.legend = FALSE) +
  labs(title = "H2O2 Stress", x = "Median Fold Change", y = "Control Noise (Normalized STD)") +
  scale_color_manual(values = c("black" = "black", "red" = "#E41A1C")) +
  theme_minimal(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16))

# --- Starvation Stress ---
ggplot(highlight_outliers(noise_STARVATION, "Starvation.FoldChange"), aes(x = Starvation.FoldChange, y = Control.STD.Normalized, label = Gene)) +
  geom_text(aes(color = color), fontface = "bold", size = 3.5, show.legend = FALSE) +
  labs(title = "Starvation", x = "Median Fold Change", y = "Control Noise (Normalized STD)") +
  scale_color_manual(values = c("black" = "black", "red" = "#E41A1C")) +
  theme_minimal(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16))


# =============================================================================
# SECTION 4: VISUALIZATION II - NOISE vs. NOISE
# =============================================================================
# These plots directly compare the expression noise (Normalized STD) of each
# protein in control conditions versus its noise level under stress.

# --- DTT Stress ---
ggplot(highlight_outliers(noise_DTT, "DTT.STD.Normalized"), aes(x = Control.STD.Normalized, y = DTT.STD.Normalized, label = Gene)) +
  geom_text(aes(color = color), fontface = "bold", size = 3.5, show.legend = FALSE) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "blue") + # y=x line
  labs(title = "DTT Stress: Noise Comparison", x = "Control Noise", y = "DTT Noise") +
  scale_color_manual(values = c("black" = "black", "red" = "#E41A1C")) +
  theme_minimal(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16))

# --- H2O2 Stress ---
ggplot(highlight_outliers(noise_H2O2, "H2O2.STD.Normalized"), aes(x = Control.STD.Normalized, y = H2O2.STD.Normalized, label = Gene)) +
  geom_text(aes(color = color), fontface = "bold", size = 3.5, show.legend = FALSE) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "blue") + # y=x line
  labs(title = "H2O2 Stress: Noise Comparison", x = "Control Noise", y = "H2O2 Noise") +
  scale_color_manual(values = c("black" = "black", "red" = "#E41A1C")) +
  theme_minimal(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16))

# --- Starvation Stress ---
ggplot(highlight_outliers(noise_STARVATION, "Starvation.STD.Normalized"), aes(x = Control.STD.Normalized, y = Starvation.STD.Normalized, label = Gene)) +
  geom_text(aes(color = color), fontface = "bold", size = 3.5, show.legend = FALSE) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "blue") + # y=x line
  labs(title = "Starvation: Noise Comparison", x = "Control Noise", y = "Starvation Noise") +
  scale_color_manual(values = c("black" = "black", "red" = "#E41A1C")) +
  theme_minimal(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16))