# =============================================================================
# Author:            Jaime Martínez Cazón
#
# Description:
# This script generates a scatter plot to visually synthesize the results from
# the candidate selection analysis. The plot maps each protein based on two key
# metrics:
#   - X-axis: The fold change in expression under nutrient starvation (log scale).
#   - Y-axis: The intrinsic expression noise under normal conditions (log scale).
#
# This visualization allows for the rapid identification of proteins in distinct
# regions of interest:
#   - Right side: Highly induced proteins.
#   - Left side: Highly repressed proteins.
#   - Top region: Proteins with high intrinsic noise.
#
# The script uses the `ggrepel` library to intelligently label key candidate
# genes in these regions, ensuring that the labels are readable and do not
# overlap, which is critical for a clear final figure.
#
# =============================================================================


# =============================================================================
# LOAD LIBRARIES
# =============================================================================
library(ggplot2)      # For the core plotting functionality.
library(dplyr)        # For data filtering and manipulation.
library(tidyverse)    # A collection of R packages for data science, used here for its general utilities.
library(ggrepel)      # For intelligent, non-overlapping text labels in ggplot.


# =============================================================================
# SECTION 1: DATA LOADING AND INTEGRATION
# =============================================================================

# --- 1.1 Load Processed Datasets ---
# Load the baseline noise and starvation stress datasets.
protein_stats_cyto <- read.csv("data/Single cell/protein_stats_cyto.csv")
noise_STARVATION <- readRDS("data/Different conditions/noise_STARVATION.rds")

# --- 1.2 Merge Datasets for Plotting ---
# Combine the two datasets using the ORF as the key to create a single,
# comprehensive data frame for visualization.
data <- merge(protein_stats_cyto, noise_STARVATION, by = "ORF", all = TRUE)


# =============================================================================
# SECTION 2: DEFINE GENE LISTS FOR TARGETED LABELING
# =============================================================================
# Manually define lists of specific genes of interest that were identified
# in the previous analysis script. These genes will be explicitly labeled on the plot.

# Genes with high induction (up-regulation) during starvation.
genes_up <- c('PDC6', 'MET14', 'VPS25','HSP12', 'TFS1', 'PBI2', 'AIM2',
              'NEW1', 'YJR096W', 'DUR1', 'MET5')

# Genes with high repression (down-regulation) during starvation.
genes_down <- c('RPS30B', 'EFT1', 'RPL40B', 'DBP2', 'LYS20', 'RPS6B',
                'ARO10', 'TEF4', 'ARO9', 'RPS25B', 'RPS24B')

# Genes with high intrinsic noise in normal conditions.
genes_top <- c('CTT1', 'YNR034W-A', 'SBP1', 'RNR4', 'HSP12', 'GRE1',
               'PGM2', 'LYS20', 'DDR48', 'EFM1', 'MSC1')
# Add a specific ORF that may be missing a gene name but is still of interest.
gene_top_orf <- c('YGR088W')


# =============================================================================
# SECTION 3: GENERATE THE SCATTER PLOT
# =============================================================================

# --- 3.1 Create the Base Scatter Plot ---
# Initialize the ggplot object.
# - The aesthetics map Fold Change to the x-axis and Normalized STD to the y-axis.
# - `geom_point()` creates the scatter plot of all proteins.
# - `scale_x_log10()` and `scale_y_log10()` transform the axes to a logarithmic
#   scale, which is essential for visualizing data that spans several orders of magnitude.
p <- ggplot(data, aes(x = Starvation.FoldChange, y = STD.MeanNormalized, label = Gene)) +
  geom_point(alpha = 0.5, color = "grey40") + # Make points slightly transparent
  scale_x_log10(
    breaks = c(0.1, 0.5, 1, 2, 5, 10), # Custom breaks for clarity
    labels = as.character(c(0.1, 0.5, 1, 2, 5, 10))
  ) +
  scale_y_log10() +
  theme_classic(base_size = 14) +
  labs(
    title = "Expression Noise vs. Starvation Response",
    x = "Fold Change under Starvation (log scale)",
    y = "Baseline Expression Noise (log scale)"
  )

# --- 3.2 Add Repulsive Labels for Different Gene Groups ---
# Layering multiple `geom_text_repel` calls allows for fine-tuned control over
# the placement of labels for each specific group of interest.

p_labeled <- p +
  # Label 1: Upregulated genes (positioned to the right).
  # `nudge_x` pushes labels horizontally, and `direction` constrains movement.
  geom_text_repel(
    data          = filter(data, Gene %in% genes_up | ORF %in% gene_top_orf),
    force         = 0.5,
    max.overlaps  = Inf,
    nudge_x       = 0.7,
    segment.size  = 0.3,
    direction     = "y",
    hjust         = 0,
    segment.color = "grey50"
  ) +
  # Label 2: Downregulated genes (positioned to the left).
  geom_text_repel(
    data          = filter(data, Gene %in% genes_down),
    force         = 0.5,
    max.overlaps  = Inf,
    nudge_x       = -0.7,
    segment.size  = 0.3,
    direction     = "y",
    hjust         = 1,
    segment.color = "grey50"
  ) +
  # Label 3: High-noise genes (positioned towards the top).
  geom_text_repel(
    data          = filter(data, Gene %in% genes_top),
    force         = 0.5,
    max.overlaps  = Inf,
    nudge_y       = 0.7,
    segment.size  = 0.3,
    direction     = "x",
    vjust         = 0,
    segment.color = "grey50"
  )

# --- 3.3 Display the Final Plot ---
print(p_labeled)