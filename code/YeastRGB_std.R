# =============================================================================
# Author:            Jaime Martínez Cazón
#
# Description:
# This script processes single-cell fluorescence microscopy data for S. cerevisiae
# from the YeastRGB database. The primary goal is to quantify protein expression
# noise for a large set of proteins identified by their Open Reading Frame (ORF).
#
# The script performs the following key steps:
# 1. Loads raw fluorescence intensity data.
# 2. Applies a quality control filter to remove unreliable measurements,
#    excluding proteins with a low cell count or a high proportion of
#    low-fluorescence cells.
# 3. For each valid protein (ORF), it calculates key statistical metrics of
#    its expression distribution, such as mean, median, standard deviation (STD),
#    and interquartile range (IQR). Noise is quantified using normalized metrics
#    (e.g., STD / mean).
# 4. Filters proteins based on functional annotations (e.g., "stress", "pH")
#    to create functionally relevant subsets.
# 5. Generates visualizations to compare the noise distribution of specific
#    functional groups against the global proteome distribution.
#
# =============================================================================


# =============================================================================
# LOAD LIBRARIES
# =============================================================================
library(ggplot2)      # For creating plots and visualizations.
library(dplyr)        # For data manipulation and transformation.
library(stringr)      # For string and pattern matching operations.
library(stringdist)   # For calculating string distances (used in an exploratory section).


# =============================================================================
# SECTION 1: DATA PROCESSING AND STATISTICAL ANALYSIS
# =============================================================================

# --- 1.1 Load and Process Raw Fluorescence Data ---

# Load the main dataset containing single-cell fluorescence measurements.
# The path is set to a relative directory for repository compatibility.
yeastRGB_data <- read.csv("data/Single cell/C_SWAT_tab_data.csv")

# Initialize a list to store fluorescence intensity values for each protein (ORF).
protein_intensities <- list()

# Iterate through each unique ORF and store its corresponding GFP intensity values.
for (orf_value in unique(yeastRGB_data$ORF)) {
  protein_intensities[[as.character(orf_value)]] <- yeastRGB_data[yeastRGB_data$ORF == orf_value, "GFP_int_b5", drop = FALSE]
}


# --- 1.2 Calculate Expression Statistics with Quality Control ---

# Initialize a list to store the calculated statistics for each protein.
protein_stats_list <- list()

# Loop through each protein to calculate statistics and apply quality filters.
for (orf_value in names(protein_intensities)) {
  intensities <- protein_intensities[[orf_value]]
  
  # Proceed only if the protein has associated measurements.
  if (nrow(intensities) > 0) {
    # Define a fluorescence threshold for quality control.
    fluorescence_threshold <- 30
    
    # Calculate the number and fraction of cells below the threshold.
    low_signal_cells <- sum(intensities$GFP_int_b5 <= fluorescence_threshold)
    total_cells <- nrow(intensities)
    high_signal_cells <- total_cells - low_signal_cells
    fraction_low_signal <- low_signal_cells / total_cells
    
    # Quality Control Filter: A protein is considered for analysis only if
    # less than 40% of its cells have low signal AND it has at least 100
    # cells with a signal above the threshold. This ensures data robustness.
    if (fraction_low_signal <= 0.4 && high_signal_cells > 100) {
      
      # Filter out low-signal cells before calculating statistics.
      intensities_filtered <- intensities[intensities$GFP_int_b5 >= fluorescence_threshold, ]
      
      # Calculate descriptive statistics for the filtered distribution.
      mean_value <- mean(intensities_filtered$GFP_int_b5, na.rm = TRUE)
      median_value <- median(intensities_filtered$GFP_int_b5, na.rm = TRUE)
      std_value <- sd(intensities_filtered$GFP_int_b5, na.rm = TRUE)
      
      # Calculate the interquartile range (IQR), a robust measure of spread.
      Q1 <- quantile(intensities_filtered$GFP_int_b5, 0.25, na.rm = TRUE)
      Q3 <- quantile(intensities_filtered$GFP_int_b5, 0.75, na.rm = TRUE)
      iqr_value <- Q3 - Q1
      
      # Calculate normalized metrics for expression noise (variability).
      # These are unitless and allow for comparison across proteins.
      normalized_std_mean <- ifelse(mean_value != 0, std_value / mean_value, NA)
      normalized_std_median <- ifelse(median_value != 0, std_value / median_value, NA)
      normalized_iqr_median <- ifelse(median_value != 0, iqr_value / median_value, NA)
      
      # Identify outliers using the 1.5*IQR rule and calculate their fraction.
      lower_bound <- Q1 - 1.5 * iqr_value
      upper_bound <- Q3 + 1.5 * iqr_value
      num_outliers <- sum(intensities_filtered$GFP_int_b5 < lower_bound | intensities_filtered$GFP_int_b5 > upper_bound, na.rm = TRUE)
      fraction_outliers <- num_outliers / nrow(intensities_filtered)
      
      # Store all calculated statistics in a data frame.
      protein_stats_list[[orf_value]] <- data.frame(
        ORF = orf_value,
        MEAN = mean_value,
        MEDIAN = median_value,
        STD = std_value,
        IQR = iqr_value,
        STD.MeanNormalized = normalized_std_mean,
        STD.MedianNormalized = normalized_std_median,
        IQR.MedianNormalized = normalized_iqr_median,
        FractionOutliers = fraction_outliers
      )
    }
  }
}

# Combine the list of data frames into a single, comprehensive data frame.
protein_stats <- do.call(rbind, protein_stats_list)

# Save the resulting statistics to a CSV file for future use.
write.csv(protein_stats, "data/Single cell/protein_stats.csv", row.names = FALSE)


# =============================================================================
# SECTION 2: FUNCTIONAL ANNOTATION AND FILTERING
# =============================================================================

# Load supplementary data containing protein names, descriptions, and keywords.
protein_info <- read.csv("data/Single cell/summarySC.csv")

# Define columns that contain functional information.
functional_columns <- c("desc", "COMMENTS", "KEYWORDS")

# Define a list of keywords to identify proteins related to specific biological processes.
keywords <- c("pH", "stress", "transporter", "DNA damage", "ATP synthesis", "TCA cycle", "Proteasome", "Vacuolar acidification")

# Initialize a list to store the subsets of proteins for each keyword.
key_list <- list()

# Loop over keywords to filter and create datasets for each functional group.
for (keyword in keywords) {
  # Filter the protein info dataframe by searching for the keyword across functional columns.
  filtered_data <- protein_info[
    grepl(keyword, protein_info$desc, ignore.case = TRUE) | 
      grepl(keyword, protein_info$COMMENTS, ignore.case = TRUE) |
      grepl(keyword, protein_info$KEYWORDS, ignore.case = TRUE),
  ]
  
  # Retain only the 'name' column, which corresponds to the ORF.
  filtered_orfs <- filtered_data["name"]
  
  # Assign the resulting list of ORFs to a dynamically named variable (e.g., "prot_pH").
  variable_name <- paste0("prot_", gsub(" ", "_", keyword)) # Replace spaces for valid names
  assign(variable_name, filtered_orfs)
  
  # Store the filtered ORF list for easy access.
  key_list[[variable_name]] <- filtered_orfs
}

# Display the names of the newly created variables containing ORF lists.
print("Created the following protein lists based on keywords:")
ls(pattern = "^prot_")


# =============================================================================
# SECTION 3: COMPARATIVE NOISE ANALYSIS AND VISUALIZATION
# =============================================================================

# Load the previously saved protein statistics file.
protein_stats <- read.csv("data/Single cell/protein_stats.csv")

# --- 3.1 Define a Protein Subset for Comparison ---

# Select a specific functional group to highlight in the plot.
# Example: Use the "pH" related proteins identified in the previous section.
# The column `name` from `prot_pH` contains the ORFs.
select_key <- prot_pH$name  

# Add a new column to the main stats dataframe to flag proteins belonging to the selected group.
protein_stats$is_key_protein <- ifelse(protein_stats$ORF %in% select_key, "Key Protein", "Other")

# --- 3.2 Generate Comparative Histogram ---

# Create a histogram to compare the noise distribution (STD.MeanNormalized) of the
# selected "Key Proteins" against all "Other" proteins.
# The plot specifically annotates HSP12 ('YFL014W'), a known stress protein.
ggplot(protein_stats, aes(x = STD.MeanNormalized, fill = is_key_protein)) +
  geom_histogram(binwidth = 0.01, aes(alpha = is_key_protein), position = "identity") +
  scale_fill_manual(values = c("Key Protein" = "#E41A1C", "Other" = "#377EB8")) + # Red/Blue palette
  scale_alpha_manual(values = c("Key Protein" = 0.8, "Other" = 0.3)) +
  labs(
    title = "Distribution of Protein Expression Noise",
    subtitle = "Highlighting pH-Related Proteins",
    x = "Expression Noise (Normalized STD)",
    y = "Number of Proteins",
    fill = "Protein Group",
    alpha = "Protein Group"
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "bottom", plot.title = element_text(face = "bold")) +
  
  # Add an arrow to highlight the position of HSP12 (ORF: YFL014W)
  geom_segment(
    data = subset(protein_stats, ORF == "YFL014W"),
    aes(x = STD.MeanNormalized, xend = STD.MeanNormalized, y = 15, yend = 1),
    arrow = arrow(length = unit(0.3, "cm")),
    color = "black",
    size = 0.7
  ) +
  
  # Add a text label for HSP12
  geom_text(
    data = subset(protein_stats, ORF == "YFL014W"),
    aes(x = STD.MeanNormalized, y = 18, label = 'HSP12'),
    stat = "identity",
    size = 4.0,
    fontface = "bold",
    color = "black"
  )