# =============================================================================
# Author:            Jaime Martínez Cazón
#
# Description:
# This script generates a multi-panel plot to visualize and compare the single-cell
# fluorescence intensity distributions for a curated list of S. cerevisiae proteins.
#
# Key features of the script include:
# 1. Data Structuring: It processes a list of data frames (`Fluo_orf`) into a
#    single long-format data frame suitable for ggplot2.
# 2. Customization: It defines specific groups of proteins, assigns them distinct
#    colors, and arranges them in a pre-determined order for clear comparison.
# 3. Advanced Plotting with ggplot2 & Patchwork: It uses `facet_wrap` to create
#    individual histograms for each protein and then employs the `patchwork`
#    library to assemble these plots into a clean, two-column layout.
#
# Pre-requisite: This script assumes that the `Fluo_orf` R object, containing
# the fluorescence data, is already loaded into the environment.
#
# =============================================================================


# =============================================================================
# LOAD LIBRARIES
# =============================================================================
library(ggplot2)      # For creating advanced data visualizations.
library(dplyr)        # For data manipulation and transformation.
library(forcats)      # For handling categorical variables (factors).
library(patchwork)    # For combining multiple ggplot plots into a single figure.


# =============================================================================
# SECTION 1: CONFIGURATION AND SETUP
# =============================================================================

# --- 1.1 Define Protein Lists, Groups, and Custom Order ---

# A comprehensive list of Open Reading Frames (ORFs) to be included in the plot.
orf_list <- c(
  "YAL049C", "YBR208C", "YDL182W", "YFL014W", "YGR088W",
  "YGR180C", "YHL034C", "YJR096W", "YJR137C", "YKL001C",
  "YLR178C", "YML128C", "YMR105C", "YPL223C", "YPL226W"
)

# Define two distinct functional or characteristic groups for coloring.
genes_orange <- c("YDL182W", "YFL014W", "YGR088W", "YGR180C", "YHL034C", "YML128C", "YMR105C", "YPL223C")
genes_purple <- c("YKL001C", "YLR178C", "YAL049C", "YPL226W", "YJR096W", "YJR137C", "YBR208C")

# Define a specific order for arranging the plots. This ensures a consistent
# and logical layout, facilitating comparison between proteins.
custom_order <- c(
  "YGR088W", "YHL034C", "YGR180C", "YFL014W", "YPL223C",
  "YMR105C", "YDL182W", "YML128C", "YKL001C", "YLR178C",
  "YAL049C", "YPL226W", "YJR096W", "YBR208C", "YJR137C"
)


# =============================================================================
# SECTION 2: DATA PREPARATION AND STRUCTURING
# =============================================================================

# --- 2.1 Validate ORFs and Create Color Mapping ---

# Ensure that we only work with ORFs that are present in the loaded `Fluo_orf` data object.
valid_orfs <- intersect(orf_list, names(Fluo_orf))

# Create a named vector that maps each valid ORF to its designated color.
# This will be used by ggplot's `scale_fill_manual`.
color_mapping <- c(
  setNames(rep("#ff7f0e", length(intersect(genes_orange, valid_orfs))), intersect(genes_orange, valid_orfs)),
  setNames(rep("#9467bd", length(intersect(genes_purple, valid_orfs))), intersect(genes_purple, valid_orfs))
)

# --- 2.2 Consolidate Data into a Single DataFrame ---

# Use `bind_rows` and `lapply` to transform the list of data frames (`Fluo_orf`)
# into a single, long-format data frame. This structure is optimal for faceting in ggplot.
df_all <- bind_rows(
  lapply(valid_orfs, function(gene) {
    dat <- Fluo_orf[[gene]]
    # Handle cases where data might be missing or empty.
    if (is.null(dat) || is.null(dat$GFP_int_b5)) return(NULL)
    vec <- dat$GFP_int_b5[!is.na(dat$GFP_int_b5)]
    if (length(vec) == 0) return(NULL)
    data.frame(ORF = gene, GFP_int_b5 = vec)
  })
)

# --- 2.3 Filter Data and Apply Custom Ordering ---

# Filter the fluorescence values to a specific range to exclude outliers and focus the visualization.
# Crucially, convert the 'ORF' column to a factor with levels defined by `custom_order`.
# This forces ggplot to render the facets in our desired sequence.
df_all <- df_all %>%
  filter(GFP_int_b5 >= 0 & GFP_int_b5 <= 2800) %>%
  mutate(ORF = factor(ORF, levels = custom_order))

# Ensure the color map only contains entries for ORFs that are actually present in the final data frame.
color_mapping <- color_mapping[names(color_mapping) %in% levels(df_all$ORF)]


# =============================================================================
# SECTION 3: PLOT GENERATION
# =============================================================================
# The plotting is done in a loop to create two separate plot objects,
# which will become the two columns in the final figure.

plot_list <- list()
num_genes_col1 <- 7
plot_genes <- list(
  custom_order[1:num_genes_col1],
  custom_order[(num_genes_col1 + 1):length(custom_order)]
)

for (i in 1:2) {
  df_sub <- df_all %>% filter(ORF %in% plot_genes[[i]])
  
  p <- ggplot(df_sub, aes(x = GFP_int_b5, fill = ORF)) +
    geom_histogram(bins = 50, position = "identity", alpha = 0.8, color = "black") +
    # Create a separate plot panel for each ORF, arranged in a single column.
    facet_wrap(~ ORF, ncol = 1, scales = "free_y") +
    # Set x-axis limits and breaks for consistency.
    scale_x_continuous(limits = c(0, 2800), breaks = seq(0, 2800, 500)) +
    # Apply the custom colors and remove the legend.
    scale_fill_manual(values = color_mapping, guide = "none") +
    theme_minimal(base_size = 14) +
    # Apply theme modifications for a clean, compact appearance.
    theme(
      axis.title.y = element_blank(),     # Remove y-axis title
      axis.text.y = element_blank(),      # Remove y-axis labels
      axis.ticks.y = element_blank(),     # Remove y-axis ticks
      axis.text.x = element_text(size = 10),
      strip.text = element_text(size = 12, face = "bold"), # Style facet titles
      panel.grid.major.x = element_line(color = "gray80", linetype = "dashed")
    )
  
  # Add the x-axis label only to the plots in the second column to avoid redundancy.
  if (i == 1) {
    p <- p + labs(x = NULL)
  } else {
    p <- p + labs(x = "GFP Intensity")
  }
  
  plot_list[[i]] <- p
}


# =============================================================================
# SECTION 4: FINAL PLOT ASSEMBLY AND DISPLAY
# =============================================================================

# Use the `patchwork` library to combine the two generated plot columns
# side-by-side into a single, cohesive figure.
final_plot <- plot_list[[1]] + plot_list[[2]]

# Display the final composite plot.
print(final_plot)