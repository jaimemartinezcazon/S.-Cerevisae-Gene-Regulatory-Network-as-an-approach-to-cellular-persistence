'''
Author: Jaime Martínez Cazón

Description:
This script visualizes the abundance of different directed triangle motifs
in the S. cerevisiae Gene Regulatory Network. It reads pre-calculated counts
for the real network and the mean/standard deviation for a surrogate model
ensemble from text files. The results are presented in a log-scale scatter
plot to compare the real network against the null model.
'''

import os
import numpy as np
import matplotlib.pyplot as plt

# =============================================================================
# SETUP: FILE PATHS AND MOTIF DEFINITIONS
# =============================================================================
DATA_DIR = "data"
REAL_DATA_FILE = os.path.join(DATA_DIR, "edge_list_stat.dat")
SURROGATE_DATA_FILE = os.path.join(DATA_DIR, "surrogate_list_stat.dat")

# Define the order and names of triangle motifs for the plot's x-axis.
# This order must match the names used in the input data files.
TRIANGLE_MOTIFS = [
    "3-cycle", "3-nocycle", "4-1biout", "4-1biin",
    "4-1biflow", "5-2bi", "6-3bi"
]

# =============================================================================
# DATA LOADING FUNCTION
# =============================================================================

def load_motif_data(real_path, surrogate_path, motifs):
    """
    Loads and parses the motif count data for both the real and surrogate networks.

    Args:
        real_path (str): Path to the real network data file.
        surrogate_path (str): Path to the surrogate model data file.
        motifs (list): An ordered list of motif names.

    Returns:
        tuple: A tuple containing (real_counts, surrogate_means, surrogate_stds).
               Returns (None, None, None) if files are not found.
    """
    if not os.path.exists(real_path) or not os.path.exists(surrogate_path):
        print("Error: One or both data files could not be found.")
        print(f"Checked for '{real_path}' and '{surrogate_path}'")
        return None, None, None

    motif_to_index = {name: i for i, name in enumerate(motifs)}

    # Initialize arrays to store data in the correct order
    real_counts = np.zeros(len(motifs))
    surrogate_means = np.zeros(len(motifs))
    surrogate_stds = np.zeros(len(motifs))

    # --- Parse real network data ---
    with open(real_path, 'r') as f:
        next(f)  # Skip header
        for line in f:
            if line.strip().startswith("total"): break
            parts = line.strip().split()
            if len(parts) >= 2:
                motif_name = parts[0]
                if motif_name in motif_to_index:
                    try:
                        count = int(parts[-1])
                        real_counts[motif_to_index[motif_name]] = count
                    except ValueError:
                        print(f"Warning: Could not parse count for motif '{motif_name}' in real data.")

    # --- Parse surrogate model data ---
    with open(surrogate_path, 'r') as f:
        next(f)  # Skip header
        for line in f:
            if line.strip().startswith("total"): break
            parts = line.strip().rsplit(None, 3) # Split from right
            if len(parts) == 4 and parts[2] == '+/-':
                motif_name = parts[0]
                if motif_name in motif_to_index:
                    try:
                        mean = float(parts[1])
                        std = float(parts[3])
                        idx = motif_to_index[motif_name]
                        surrogate_means[idx] = mean
                        surrogate_stds[idx] = std
                    except ValueError:
                        print(f"Warning: Could not parse mean/std for motif '{motif_name}' in surrogate data.")

    return real_counts, surrogate_means, surrogate_stds

# =============================================================================
# PLOTTING FUNCTION
# =============================================================================

def create_motif_plot(real_counts, surrogate_means, surrogate_stds, motifs):
    """
    Generates and displays the scatter plot for triangle motif comparison.

    Args:
        real_counts (np.ndarray): Counts for each motif in the real network.
        surrogate_means (np.ndarray): Mean counts for each motif in the null model.
        surrogate_stds (np.ndarray): Std dev for each motif in the null model.
        motifs (list): Ordered list of motif names.
    """
    plt.figure(figsize=(12, 8))

    x_positions = np.arange(len(motifs))
    x_tick_labels = [str(i + 1) for i in range(len(motifs))]

    # --- Plot Real Network Data ---
    # We only plot points with counts > 0 for clarity on a log scale
    real_indices_to_plot = np.where(real_counts > 0)[0]
    if len(real_indices_to_plot) > 0:
        plt.scatter(x_positions[real_indices_to_plot], real_counts[real_indices_to_plot],
                    marker='s', color='#bcbd22', s=150, zorder=10,
                    label='Real Network')

    # --- Plot Surrogate Model Data with Error Bars ---
    # We only plot points with mean counts > 0
    surrogate_indices_to_plot = np.where(surrogate_means > 0)[0]
    if len(surrogate_indices_to_plot) > 0:
        plt.errorbar(x_positions[surrogate_indices_to_plot], surrogate_means[surrogate_indices_to_plot],
                     yerr=surrogate_stds[surrogate_indices_to_plot],
                     fmt='o', color='black', markersize=10, capsize=5,
                     alpha=0.8, linestyle='None', label='Configuration Model')

    # --- Plot Formatting ---
    plt.yscale('log')
    plt.title("Triangle Motif Abundance", fontsize=28, fontweight='bold')
    plt.xlabel("Triangle Motif Category", fontsize=22)
    plt.ylabel("Number of Triangles (#)", fontsize=22)

    plt.xticks(x_positions, x_tick_labels, fontsize=18)
    plt.yticks(fontsize=18)

    plt.legend(fontsize=16)
    plt.tight_layout()
    plt.show()

# =============================================================================
# SCRIPT EXECUTION
# =============================================================================

if __name__ == "__main__":
    # 1. Load the pre-calculated data
    real_data, surrogate_mean_data, surrogate_std_data = load_motif_data(
        REAL_DATA_FILE, SURROGATE_DATA_FILE, TRIANGLE_MOTIFS
    )

    # 2. Check if data was loaded successfully and create the plot
    if real_data is not None:
        print("Data loaded successfully. Generating plot...")
        create_motif_plot(real_data, surrogate_mean_data, surrogate_std_data, TRIANGLE_MOTIFS)
    else:
        print("Could not generate plot due to missing data.")