'''
Author: Jaime Martínez Cazón

Description:
This script performs a detailed analysis of the degree distributions (in, out,
and total) of the S. cerevisiae Gene Regulatory Network. It visualizes the
Complementary Cumulative Distribution Function (CCDF) for each degree type on a
log-log scale. The script applies exponential binning to the raw CCDF data to
reveal underlying trends and performs a power-law fit on the total degree
distribution to estimate the scaling exponent gamma, a key characteristic of
scale-free networks.
'''

import os
import numpy as np
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
from scipy.stats import linregress

# =============================================================================
# SETUP: FILE PATHS
# =============================================================================
DATA_DIR = "data"
EDGE_LIST_PATH = os.path.join(DATA_DIR, "giantC_edge_list.csv")

# =============================================================================
# UTILITY FUNCTIONS
# =============================================================================

def calculate_ccdf(degrees):
    """
    Calculates the CCDF points (k, P(K>=k)) for a list of degrees.

    Args:
        degrees (list or np.ndarray): A list of degree values.

    Returns:
        tuple: A tuple of numpy arrays (unique_degrees, ccdf_values).
    """
    if not isinstance(degrees, np.ndarray):
        degrees = np.array(degrees)
    
    unique_degrees, counts = np.unique(degrees, return_counts=True)
    cumulative_counts = np.cumsum(counts[::-1])[::-1]
    ccdf = cumulative_counts / len(degrees)
    
    return unique_degrees, ccdf

def exponential_binning(k_data, val_data, base=1.5):
    """
    Applies exponential binning to (k, value) data points.

    Args:
        k_data (np.ndarray): Array of degrees.
        val_data (np.ndarray): Array of corresponding values (e.g., CCDF).
        base (float): The base for the exponential bins.

    Returns:
        tuple: (binned_k_mean, binned_val_mean, binned_val_std)
    """
    if len(k_data) == 0:
        return np.array([]), np.array([]), np.array([])
        
    min_k, max_k = np.min(k_data), np.max(k_data)
    
    # Define bin edges
    current_edge = float(min_k)
    bin_edges = [current_edge]
    while current_edge <= max_k:
        current_edge *= base
        bin_edges.append(current_edge)

    # Bin data and calculate statistics
    binned_k, binned_val, binned_std = [], [], []
    for i in range(len(bin_edges) - 1):
        low, high = bin_edges[i], bin_edges[i+1]
        indices = np.where((k_data >= low) & (k_data < high))[0]
        
        if len(indices) > 0:
            binned_k.append(np.mean(k_data[indices]))
            binned_val.append(np.mean(val_data[indices]))
            binned_std.append(np.std(val_data[indices]))
            
    return np.array(binned_k), np.array(binned_val), np.array(binned_std)


# =============================================================================
# PLOTTING FUNCTION
# =============================================================================

def plot_ccdf_distribution(ax, degrees, title, fit_power_law=False):
    """
    Plots the CCDF, binned data, and optional power-law fit for a degree distribution.

    Args:
        ax (matplotlib.axes.Axes): The matplotlib axes object to plot on.
        degrees (list): The list of degree values for the distribution.
        title (str): The title for the plot.
        fit_power_law (bool): Whether to perform and plot a power-law fit.
    """
    # --- 1. Calculate CCDF and filter for plotting ---
    k_unique, ccdf_values = calculate_ccdf(degrees)
    valid_indices = (k_unique > 0) & (ccdf_values > 0)
    k_plot, ccdf_plot = k_unique[valid_indices], ccdf_values[valid_indices]
    
    if len(k_plot) == 0:
        ax.text(0.5, 0.5, "No data to plot", ha='center', va='center', transform=ax.transAxes)
        return

    # --- 2. Plot raw CCDF data ---
    ax.scatter(k_plot, ccdf_plot, marker='o', s=30, alpha=0.2, color='#8F9058', label='Raw CCDF Data', zorder=1)
    
    # --- 3. Exponential Binning and Plotting ---
    k_binned, ccdf_binned, ccdf_std = exponential_binning(k_plot, ccdf_plot)
    ax.errorbar(k_binned, ccdf_binned, yerr=ccdf_std, fmt='s', color='#bcbd22',
                markersize=8, capsize=5, label='Binned CCDF (Mean ± SD)', zorder=10)
                
    # --- 4. Power-Law Fit (if requested) ---
    if fit_power_law and len(k_plot) > 1:
        log_k = np.log(k_plot)
        log_ccdf = np.log(ccdf_plot)
        
        try:
            slope, intercept, r_value, _, _ = linregress(log_k, log_ccdf)
            gamma = -slope + 1 # Exponent of the PDF, P(k) ~ k^-gamma
            
            fit_k_range = np.logspace(np.log10(min(k_plot)), np.log10(max(k_plot)), 100)
            fit_ccdf = np.exp(intercept) * (fit_k_range ** slope)
            
            ax.plot(fit_k_range, fit_ccdf, color='#d62728', linewidth=3,
                    label=f'Power-Law Fit (γ ≈ {gamma:.2f})', zorder=9)
            
            print(f"\n--- Power-Law Fit Results for {title} ---")
            print(f"Slope of CCDF (log-log scale): {slope:.2f}")
            print(f"Estimated PDF exponent (gamma): {gamma:.2f}")
            print(f"R-squared of log-log fit: {r_value**2:.2f}")
            print("-" * 40)
            
        except ValueError as e:
            print(f"Could not perform power-law fit for {title}: {e}")

    # --- 5. Final Plot Formatting ---
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel("Degree (k)", fontsize=18)
    ax.set_ylabel("CCDF P(K ≥ k)", fontsize=18)
    ax.set_title(title, fontsize=22, fontweight='bold')
    ax.tick_params(labelsize=16)
    ax.legend(fontsize=14)
    ax.grid(False)


# =============================================================================
# MAIN EXECUTION
# =============================================================================

if __name__ == "__main__":
    # --- 1. Load and prepare network data ---
    print("Loading network and calculating degrees...")
    edge_list = pd.read_csv(EDGE_LIST_PATH)
    G_full = nx.from_pandas_edgelist(
        edge_list, source='Node 1', target='Node 2', create_using=nx.DiGraph()
    )
    
    # Extract the giant component to ensure a connected graph for analysis
    giant_component_nodes = max(nx.weakly_connected_components(G_full), key=len)
    G = G_full.subgraph(giant_component_nodes).copy()
    
    # Get degree lists
    in_degrees = [d for _, d in G.in_degree()]
    out_degrees = [d for _, d in G.out_degree()]
    total_degrees = [d for _, d in G.to_undirected().degree()]

    # --- 2. Create plots ---
    print("Generating CCDF plots...")
    fig, axes = plt.subplots(1, 3, figsize=(24, 7))
    
    plot_ccdf_distribution(axes[0], in_degrees, "In-Degree CCDF")
    plot_ccdf_distribution(axes[1], out_degrees, "Out-Degree CCDF")
    plot_ccdf_distribution(axes[2], total_degrees, "Total Degree CCDF", fit_power_law=True)
    
    plt.tight_layout()
    plt.show()

    print("\nAnalysis complete.")