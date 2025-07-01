'''
Author: Jaime Martínez Cazón

Description:
This script analyzes the relationship between the clustering coefficient C(k) 
and node degree (k) for the S. cerevisiae Gene Regulatory Network. It compares 
the properties of the real network against an ensemble of surrogate networks 
(null models). The analysis includes:
1.  Calculation of the global clustering coefficient and its statistical 
    significance.
2.  Plotting C(k) vs. k using exponential binning for the real and null models.
3.  Performing a power-law fit (C(k) ~ k^gamma) to the binned real network data
    to characterize its hierarchical structure.
'''

import os
import numpy as np
import pandas as pd
import networkx as nx
from tqdm import tqdm
from collections import defaultdict
from scipy.stats import linregress
import matplotlib.pyplot as plt

# =============================================================================
# SETUP: DIRECTORIES AND FILE PATHS
# =============================================================================
DATA_DIR = "data"
EDGE_LIST_PATH = os.path.join(DATA_DIR, "giantC_edge_list.csv")
NULL_MODELS_DIR = os.path.join(DATA_DIR, "Null Models")
N_NULL_MODELS = 1000  # Number of surrogate models to analyze

# =============================================================================
# UTILITY FUNCTIONS
# =============================================================================

def calculate_empirical_p_value(real_value, null_distribution, direction='greater'):
    """
    Calculates the empirical p-value based on a null distribution.
    
    Args:
        real_value (float): The metric's value in the real network.
        null_distribution (list): The metric's values from the null models.
        direction (str): 'greater' if higher values are more significant, 
                         'less' otherwise.

    Returns:
        float: The empirical p-value.
    """
    n_simulations = len(null_distribution)
    if n_simulations == 0:
        return np.nan
        
    null_array = np.array(null_distribution)
    
    if direction == 'greater':
        count = np.sum(null_array >= real_value)
    else:  # direction == 'less'
        count = np.sum(null_array <= real_value)
        
    p_value = (count + 1) / (n_simulations + 1)
    return p_value


def exponential_binning(data_points, base=1.5):
    """
    Applies exponential binning to a list of (degree, value) data points.

    Args:
        data_points (list): A list of tuples (k, C(k)).
        base (float): The base for the exponential bins (e.g., 1.5).

    Returns:
        tuple: (binned_k_mean, binned_c_mean, binned_c_std)
    """
    if not data_points:
        return np.array([]), np.array([]), np.array([])

    degrees = np.array([k for k, _ in data_points])
    c_values = np.array([c for _, c in data_points])

    max_degree = np.max(degrees)
    min_degree = np.min(degrees)
    
    # Define bin edges
    current_edge = float(min_degree)
    bin_edges = [current_edge]
    while current_edge <= max_degree:
        current_edge *= base
        bin_edges.append(current_edge)

    # Bin the data and calculate statistics
    binned_k_means, binned_c_means, binned_c_stds = [], [], []
    for i in range(len(bin_edges) - 1):
        low_bound, high_bound = bin_edges[i], bin_edges[i+1]
        
        indices = np.where((degrees >= low_bound) & (degrees < high_bound))[0]
        if len(indices) > 0:
            binned_k_means.append(np.mean(degrees[indices]))
            binned_c_means.append(np.mean(c_values[indices]))
            binned_c_stds.append(np.std(c_values[indices]))

    return np.array(binned_k_means), np.array(binned_c_means), np.array(binned_c_stds)

# =============================================================================
# DATA LOADING AND PREPARATION
# =============================================================================

def load_real_network():
    """Loads the real network and extracts its giant component."""
    print("Loading real network...")
    edge_list = pd.read_csv(EDGE_LIST_PATH)
    G_full = nx.from_pandas_edgelist(edge_list, source='Node 1', target='Node 2', create_using=nx.DiGraph())
    
    # Extract the giant component from the undirected version of the graph
    undirected_components = (G_full.subgraph(c) for c in nx.connected_components(G_full.to_undirected()))
    G_real = max(undirected_components, key=len)
    G_real = G_full.subgraph(G_real.nodes()).copy() # Ensure it's a DiGraph
    
    print(f"Real network (giant component) loaded: {G_real.number_of_nodes()} nodes, {G_real.number_of_edges()} edges.")
    return G_real

def analyze_null_models():
    """Loads and analyzes the ensemble of null models."""
    print(f"\nLoading and analyzing up to {N_NULL_MODELS} null models...")
    null_global_clustering = []
    # Store all C(k) values for each degree across all models for robust binning
    null_ck_raw_data = {"all": defaultdict(list), "in": defaultdict(list), "out": defaultdict(list)}
    
    for i in tqdm(range(N_NULL_MODELS)):
        file_path = os.path.join(NULL_MODELS_DIR, f"null_model_{str(i).zfill(4)}.graphml")
        if not os.path.exists(file_path):
            continue
            
        try:
            G_null = nx.read_graphml(file_path, node_type=int)
            
            local_c_null = nx.clustering(G_null.to_undirected())
            null_global_clustering.append(np.mean(list(local_c_null.values())))

            # Store C(k) values for each degree type
            for node, c_val in local_c_null.items():
                if c_val > 0:
                    null_ck_raw_data["all"][G_null.degree(node)].append(c_val)
                    null_ck_raw_data["in"][G_null.in_degree(node)].append(c_val)
                    null_ck_raw_data["out"][G_null.out_degree(node)].append(c_val)
        except Exception as e:
            tqdm.write(f"Warning: Could not process {os.path.basename(file_path)}. Error: {e}")

    return null_global_clustering, null_ck_raw_data


# =============================================================================
# PLOTTING AND ANALYSIS
# =============================================================================

def plot_ck_distribution(real_data, null_data, degree_type):
    """
    Creates a C(k) vs. k plot with exponential binning and power-law fit.
    
    Args:
        real_data (dict): Data for the real network.
        null_data (defaultdict): Raw C(k) data for the null models.
        degree_type (str): Type of degree to plot ('all', 'in', 'out').
    """
    plt.figure(figsize=(9, 7))
    ax = plt.gca()

    # --- 1. Process and plot real network data ---
    real_raw_points = [(k, c) for k, c in zip(real_data[f'degree_{degree_type}'], real_data['local_clustering']) if k > 0 and c > 0]
    if not real_raw_points:
        print(f"Warning: No valid data points for real network ({degree_type} degree).")
        return

    # Exponential binning for real data
    k_real_binned, c_real_binned, c_real_std = exponential_binning(real_raw_points)
    
    # Plot raw and binned data
    ax.scatter([k for k,c in real_raw_points], [c for k,c in real_raw_points],
               color='#8F9058', marker='o', s=20, alpha=0.1, label='Real Network (Raw)', zorder=1)
    ax.errorbar(k_real_binned, c_real_binned, yerr=c_real_std, fmt='s', color='#bcbd22',
                markersize=8, label='Real Network (Binned)', zorder=10)

    # --- 2. Process and plot null model data ---
    null_raw_points = [(k, c) for k, c_list in null_data[degree_type].items() for c in c_list if k > 0]
    if null_raw_points:
        k_null_binned, c_null_binned, c_null_std = exponential_binning(null_raw_points)
        ax.errorbar(k_null_binned, c_null_binned, yerr=c_null_std, fmt='o', color='black',
                    markersize=5, alpha=0.8, label='Null Model (Binned)', zorder=5)

    # --- 3. Power-law fit for binned real data ---
    if degree_type != "in" and len(k_real_binned) > 1:
        valid_indices = c_real_binned > 0
        if np.sum(valid_indices) > 1:
            log_k = np.log(k_real_binned[valid_indices])
            log_c = np.log(c_real_binned[valid_indices])
            slope, intercept, _, _, _ = linregress(log_k, log_c)
            
            fit_k = np.linspace(min(k_real_binned), max(k_real_binned), 100)
            fit_c = np.exp(intercept) * (fit_k ** slope)
            ax.plot(fit_k, fit_c, color='#d62728', linewidth=3,
                    label=fr'Fit ($\gamma={slope:.2f}$)', zorder=9)

    # --- 4. Final plot formatting ---
    ax.set_xscale('log')
    if degree_type != "in":
        ax.set_yscale('log')

    ax.set_xlabel("Degree (k)", fontsize=18)
    ax.set_ylabel("Clustering C(k)", fontsize=18)
    ax.set_title(f"C(k) vs. {degree_type.capitalize()} Degree", fontsize=22, fontweight='bold')
    ax.tick_params(labelsize=16)
    ax.legend(fontsize=14)
    plt.tight_layout()
    plt.show()

# =============================================================================
# MAIN EXECUTION
# =============================================================================

if __name__ == "__main__":
    # --- 1. Load and analyze real network ---
    G_real = load_real_network()
    
    # Calculate properties for the real network
    real_network_data = {
        'local_clustering': list(nx.clustering(G_real.to_undirected()).values()),
        'global_clustering': np.mean(list(nx.clustering(G_real.to_undirected()).values())),
        'degree_all': [d for _, d in G_real.degree()],
        'degree_in': [d for _, d in G_real.in_degree()],
        'degree_out': [d for _, d in G_real.out_degree()]
    }
    
    # --- 2. Load and analyze null models ---
    null_global_clustering, null_ck_raw_data = analyze_null_models()
    
    # --- 3. Perform statistical significance test for global clustering ---
    print("\n--- Global Clustering Significance Analysis ---")
    if null_global_clustering:
        p_value = calculate_empirical_p_value(
            real_network_data['global_clustering'], 
            null_global_clustering
        )
        print(f"Real Network Global Clustering: {real_network_data['global_clustering']:.4f}")
        print(f"Null Model Mean Clustering:     {np.mean(null_global_clustering):.4f} ± {np.std(null_global_clustering):.4f}")
        print(f"Empirical P-value:              {p_value:.4f}")
        if p_value < 0.05:
            print("Result is statistically significant (p < 0.05).")
        else:
            print("Result is not statistically significant (p >= 0.05).")
    else:
        print("No null model data available for significance testing.")
    print("-" * 50)
    
    # --- 4. Generate plots for each degree type ---
    print("\nGenerating C(k) vs. Degree plots...")
    for deg_type in ["all", "in", "out"]:
        plot_ck_distribution(real_network_data, null_ck_raw_data, deg_type)