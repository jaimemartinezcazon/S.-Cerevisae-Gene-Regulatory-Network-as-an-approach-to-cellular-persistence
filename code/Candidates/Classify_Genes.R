# =============================================================================
# Author:            Jaime Martínez Cazón
#
# Description:
# This script integrates data from two distinct high-throughput microscopy
# datasets (YeastRGB and LoQate) to identify candidate proteins in S. cerevisiae
# that exhibit extreme expression dynamics. The goal is to nominate proteins
# for further study based on two primary criteria:
#
# 1. High Intrinsic Noise: Identifies proteins that display the highest
#    expression variability (noise) under normal, unstressed conditions.
#    These are the top 1% of proteins based on their normalized standard deviation.
#
# 2. Extreme Response to Stress: Identifies proteins whose expression levels
#    change most dramatically in response to nutrient starvation. This includes
#    both the most highly induced (top 1%) and the most highly repressed (top 1%)
#    proteins.
#
# By merging these datasets, the script creates a unified resource for ranking
# and selecting proteins based on these quantitative criteria.
#
# =============================================================================


# =============================================================================
# LOAD LIBRARIES
# =============================================================================
library(ggplot2)      # For data visualization.
library(dplyr)        # For data manipulation and transformation.
library(stringr)      # For string and pattern matching operations.


# =============================================================================
# SECTION 1: DATA LOADING AND INTEGRATION
# =============================================================================

# --- 1.1 Load Processed Datasets ---
# Load the datasets containing baseline noise statistics, stress response data,
# and a mapping file for ORF-to-gene names. Paths are relative.

# Baseline noise data for cytoplasmic proteins (from YeastRGB).
protein_stats_cyto <- read.csv("data/Single cell/protein_stats_cyto.csv")

# ORF-to-Gene name mapping file.
names <- readRDS("data/Single cell/ORF_genes.rds")

# Starvation stress response data (from LoQate).
noise_STARVATION <- readRDS("data/Different conditions/noise_STARVATION.rds")


# --- 1.2 Merge Datasets for Integrated Analysis ---

# First, add gene names to the baseline noise statistics data.
protein_stats_cyto <- merge(names, protein_stats_cyto, by = "ORF", all = FALSE)

# Now, merge the baseline data with the starvation stress data using the ORF
# as the common key. An outer join (all = TRUE) is used to retain all proteins
# from both datasets, even if one is missing from the other.
integrated_data <- merge(protein_stats_cyto, noise_STARVATION, by = "ORF", all = TRUE)


# =============================================================================
# SECTION 2: CANDIDATE SELECTION I - HIGH BASELINE NOISE
# =============================================================================
# This section identifies proteins with the highest intrinsic expression noise
# under normal growth conditions. We use the 99th percentile as the threshold
# to define the "high-noise" category.

# --- 2.1 Define High-Noise Thresholds ---
# Calculate the 99th percentile for noise, using both mean-normalized and
# median-normalized standard deviation for robustness.
threshold_noise_mean <- quantile(protein_stats_cyto$STD.MeanNormalized, probs = 0.99, na.rm = TRUE)
threshold_noise_median <- quantile(protein_stats_cyto$STD.MedianNormalized, probs = 0.99, na.rm = TRUE)

# --- 2.2 Filter for High-Noise Candidates ---
# Select proteins whose noise value is greater than or equal to the threshold.
top_noise_mean <- protein_stats_cyto[
  which(protein_stats_cyto$STD.MeanNormalized >= threshold_noise_mean),
  c("ORF", "Gene", "STD.MeanNormalized")
]

top_noise_median <- protein_stats_cyto[
  which(protein_stats_cyto$STD.MedianNormalized >= threshold_noise_median),
  c("ORF", "Gene", "STD.MedianNormalized")
]

# --- 2.3 Rank Candidates ---
# Order the candidate lists from highest to lowest noise.
top_noise_mean_ranked <- top_noise_mean[order(top_noise_mean$STD.MeanNormalized, decreasing = TRUE), ]
top_noise_median_ranked <- top_noise_median[order(top_noise_median$STD.MedianNormalized, decreasing = TRUE), ]


# =============================================================================
# SECTION 3: CANDIDATE SELECTION II - EXTREME FOLD CHANGE UNDER STARVATION
# =============================================================================
# This section identifies proteins that show the most significant changes in
# expression levels (up or down) during nutrient starvation.

# --- 3.1 Candidates with Highest Induction (Upregulation) ---

# Define the threshold for high induction as the 99th percentile of fold change.
threshold_induced <- quantile(
  integrated_data$Starvation.FoldChange,
  probs = 0.99,
  na.rm = TRUE
)

# Filter for proteins with a fold change greater than or equal to this threshold.
top_induced_candidates <- integrated_data[
  which(integrated_data$Starvation.FoldChange >= threshold_induced),
  c("ORF", "Gene.y", "Starvation.FoldChange") # Gene.y corresponds to the gene name from the starvation dataset
]

# Rank the candidates from most induced to least.
top_induced_ranked <- top_induced_candidates[order(top_induced_candidates$Starvation.FoldChange, decreasing = TRUE), ]


# --- 3.2 Candidates with Highest Repression (Downregulation) ---

# Define the threshold for high repression as the 1st percentile of fold change.
threshold_repressed <- quantile(
  integrated_data$Starvation.FoldChange,
  probs = 0.01,
  na.rm = TRUE
)

# Filter for proteins with a fold change less than or equal to this threshold.
top_repressed_candidates <- integrated_data[
  which(integrated_data$Starvation.FoldChange <= threshold_repressed),
  c("ORF", "Gene.y", "Starvation.FoldChange")
]

# Rank the candidates from most repressed to least (ascending order).
top_repressed_ranked <- top_repressed_candidates[order(top_repressed_candidates$Starvation.FoldChange, decreasing = FALSE), ]


# =============================================================================
# SECTION 4: DISPLAY RESULTS
# =============================================================================
cat("--- Top High-Noise Candidates (Mean-Normalized) ---\n")
print(top_noise_mean_ranked)

cat("\n--- Top High-Noise Candidates (Median-Normalized) ---\n")
print(top_noise_median_ranked)

cat("\n--- Top Induced Candidates under Starvation ---\n")
print(top_induced_ranked)

cat("\n--- Top Repressed Candidates under Starvation ---\n")
print(top_repressed_ranked)