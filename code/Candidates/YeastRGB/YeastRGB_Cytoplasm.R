# =============================================================================
# Author:            Jaime Martínez Cazón
#
# Description:
# This script is designed to filter and analyze a dataset of protein expression
# noise in S. cerevisiae. The primary goal is to isolate proteins with specific
# subcellular localizations and functional roles, and then to compare their
# expression noise profiles against the broader proteome.
#
# The script follows these main steps:
# 1. Data Loading: It loads two datasets: one with pre-calculated protein
#    expression statistics (mean, noise, etc.) and another with detailed
#    protein annotations (e.g., subcellular location, function).
# 2. Subcellular Location Filtering: It first filters the dataset to retain
#    only proteins localized to the cytoplasm/cytosol. This step is crucial
#    for data quality, as 2D fluorescence microscopy measurements are most
#    accurate for proteins that are homogeneously distributed within the cell volume.
# 3. Functional Keyword Filtering: Using the cytoplasm-filtered list, it
#    further isolates subsets of proteins based on functional keywords
#    (e.g., "stress", "DNA damage").
# 4. Comparative Visualization: It generates a histogram to visually compare
#    the expression noise distribution of a selected functional group against
#    all other cytoplasmic proteins, highlighting key differences.
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

# Load the protein expression statistics and annotation files.
# Paths are set to relative directories for repository compatibility.
protein_stats <- read.csv("data/Single cell/protein_stats.csv")
protein_info <- read.csv("data/Single cell/summarySC.csv")


# =============================================================================
# SECTION 2: FILTERING BY SUBCELLULAR LOCATION
# =============================================================================
# This section filters proteins to retain those primarily located in the
# cytoplasm. This is a critical quality control step because the 2D microscopy
# technique used assumes a homogeneous protein distribution within the cell,
# which is most applicable to cytosolic proteins. Proteins with unknown
# location (NA) are also retained for initial consideration.

# Create a subset of protein annotations for cytoplasmic/cytosolic proteins.
cytoplasm_list <- protein_info[
  grepl("cytoplasm", protein_info$SUBCELLULAR.LOCATIONS, ignore.case = TRUE) |
    grepl("cytosol", protein_info$SUBCELLULAR.LOCATIONS, ignore.case = TRUE) |
    is.na(protein_info$SUBCELLULAR.LOCATIONS),
]

# Extract the unique ORF names from the filtered list.
cytoplasm_orfs <- unique(cytoplasm_list$name)


# --- (Optional) Save the filtered statistics dataset ---
# The following block creates and saves a new CSV file containing statistics
# for cytoplasmic proteins only. It is commented out but can be run to
# generate this useful subset for other analyses.
'
# Filter the main statistics dataframe to keep only cytoplasmic proteins.
protein_stats_cyto <- protein_stats[protein_stats$ORF %in% cytoplasm_orfs, ]

# Save the filtered data.
write.csv(protein_stats_cyto, "data/Single cell/protein_stats_cyto.csv", row.names = FALSE)
'


# =============================================================================
# SECTION 3: FILTERING BY FUNCTIONAL KEYWORDS
# =============================================================================
# This section uses the cytoplasm-filtered list to create smaller, functionally
# coherent groups of proteins based on keywords related to biological processes.

# Define columns in the annotation file that contain functional information.
functional_columns <- c("desc", "COMMENTS", "KEYWORDS")

# Define a list of keywords for filtering.
keywords <- c("pH", "stress", "transporter", "DNA damage", "ATP synthesis", "TCA cycle", "Proteasome", "Vacuolar acidification")

# Initialize a list to store the resulting ORF lists for each keyword.
key_list <- list()

# Loop through each keyword to create and store a filtered dataset.
for (keyword in keywords) {
  # Filter the cytoplasmic protein list by searching for the keyword.
  filtered_data <- cytoplasm_list[
    grepl(keyword, cytoplasm_list$desc, ignore.case = TRUE) |
      grepl(keyword, cytoplasm_list$COMMENTS, ignore.case = TRUE) |
      grepl(keyword, cytoplasm_list$KEYWORDS, ignore.case = TRUE),
  ]
  
  # Retain only the 'name' column (ORFs).
  filtered_orfs <- filtered_data["name"]
  
  # Create a clean variable name (e.g., "prot_DNA_damage").
  variable_name <- paste0("prot_", gsub(" ", "_", keyword))
  
  # Assign the filtered ORF list to the new variable.
  assign(variable_name, filtered_orfs)
  
  # Store the list for easy access.
  key_list[[variable_name]] <- filtered_orfs
}

# Display the names of the newly created protein group variables.
print("Created the following protein lists based on keywords:")
ls(pattern = "^prot_")


# =============================================================================
# SECTION 4: COMPARATIVE NOISE VISUALIZATION
# =============================================================================

# --- 4.1 Define Protein Subset and Prepare Data for Plotting ---

# Select a specific functional group to highlight in the plot.
# Here, we use the "DNA damage" response proteins.
select_key <- prot_DNA_damage$name

# Add a new column to the main stats dataframe to flag proteins belonging
# to the selected functional group.
protein_stats$is_key_protein <- ifelse(protein_stats$ORF %in% select_key, "Key Protein", "Other")

# --- 4.2 Generate Comparative Histogram ---
# The plot displays two overlaid histograms: one for the "Key Protein" group
# and one for all "Other" proteins. This allows for a direct visual comparison
# of their expression noise distributions. The plot also annotates the
# position of HSP12, a well-known high-noise protein.

ggplot(protein_stats, aes(x = normalized_std, fill = is_key_protein)) +
  # Use position="identity" to overlay histograms, with alpha for transparency.
  geom_histogram(binwidth = 0.01, aes(alpha = is_key_protein), position = "identity") +
  # Define custom colors and transparency levels for clarity.
  scale_fill_manual(values = c("Key Protein" = "#E41A1C", "Other" = "#377EB8")) +
  scale_alpha_manual(values = c("Key Protein" = 0.8, "Other" = 0.3)) +
  labs(
    title = "Distribution of Protein Expression Noise",
    subtitle = "Highlighting DNA Damage-Related Proteins",
    x = "Expression Noise (Normalized STD)",
    y = "Number of Proteins",
    fill = "Protein Group",
    alpha = "Protein Group"
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "top") +
  
  # --- Annotation for HSP12 (YFL014W) ---
  # Add an arrow pointing to the noise value of HSP12 on the x-axis.
  geom_segment(
    data = subset(protein_stats, ORF == "YFL014W"),
    aes(x = normalized_std, xend = normalized_std, y = 15, yend = 1),
    arrow = arrow(length = unit(0.3, "cm")),
    color = "black",
    size = 0.7,
    inherit.aes = FALSE # Prevent aesthetics from being inherited from the main ggplot call.
  ) +
  # Add a text label for HSP12 above the arrow.
  geom_text(
    data = subset(protein_stats, ORF == "YFL014W"),
    aes(x = normalized_std, y = 18, label = 'HSP12'),
    stat = "identity",
    size = 4.0,
    fontface = "bold",
    color = "black",
    inherit.aes = FALSE
  )