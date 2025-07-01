'''
Author: Jaime Martínez Cazón

Description:
This script analyzes the assortativity of the S. cerevisiae Gene Regulatory
Network by calculating the Average Nearest-Neighbor degree (ANN) as a function
of node degree (k). Assortativity measures the tendency of nodes to connect
to other nodes with similar or dissimilar degrees. The script compares the
ANN(k) trend of the real network against an ensemble of surrogate networks.

Five cases are analyzed:
1. Undirected: ANN(k) vs. total degree k.
2. Directed (4 types): ANN(source_degree) vs. source_degree, where the
   neighbor's degree type can be 'in' or 'out'.
'''

import os
import pandas as pd
import numpy as np
import networkx as nx
from tqdm import tqdm
from collections import defaultdict
import matplotlib.pyplot as plt

# =============================================================================
# SETUP: FILE PATHS AND PARAMETERS
# =============================================================================
DATA_DIR = "data"
EDGE_LIST_PATH = os.path.join(DATA_DIR, "giantC_edge_list.csv")
NULL_MODELS_DIR = os.path.join(DATA_DIR, "Null Models")
N_NULL_MODELS = 1000  # Number of surrogate models to analyze

# =============================================================================
# UTILITY AND ANALYSIS FUNCTIONS
# =============================================================================

def exponential_binning(data_points, base=1.5):
    """
    Applies exponential binning to (degree, value) data points.
    (Reusing the same robust function from previous scripts).
    """
    if not data_points:
        return np.array([]), np.array([]), np.array([])
    
    k_data = np.array([k for k, _ in data_points])
    val_data = np.array([v for _, v in data_points])
    
    min_k, max_k = np.min(k_data), np.max(k_data)
    
    current_edge = float(min_k)
    bin_edges = [current_edge]
    while current_edge <= max_k:
        current_edge *= base
        bin_edges.append(current_edge)
    
    binned_k, binned_val, binned_std = [], [], []
    for i in range(len(bin_edges) - 1):
        low, high = bin_edges[i], bin_edges[i+1]
        indices = np.where((k_data >= low) & (k_data < high))[0]
        if len(indices) > 0:
            binned_k.append(np.mean(k_data[indices]))
            binned_val.append(np.mean(val_data[indices]))
            binned_std.append(np.std(val_data[indices]))
            
    return np.array(binned_k), np.array(binned_val), np.array(binned_std)


def load_real_network():
    """Loads the real network and extracts its giant component."""
    print("Loading real network...")
    edge_list = pd.read_csv(EDGE_LIST_PATH)
    G_full = nx.from_pandas_edgelist(
        edge_list, source='Node 1', target='Node 2', create_using=nx.DiGraph()
    )
    # Ensure nodes are integers
    G_full = nx.relabel_nodes(G_full, {n: int(n) for n in G_full.nodes()})
    
    giant_component_nodes = max(nx.weakly_connected_components(G_full), key=len)
    G_real = G_full.subgraph(giant_component_nodes).copy()
    
    print(f"Real network (giant component) loaded: {G_real.number_of_nodes()} nodes.")
    return G_real


def analyze_null_models_for_ann(null_dir, num_models):
    """
    Loads and analyzes the ensemble of null models to aggregate ANN data.

    Returns:
        dict: A dictionary where keys are assortativity types and values are
              lists of (degree, ANN) tuples aggregated from all models.
    """
    print(f"\nAnalyzing up to {num_models} null models for ANN...")
    null_ann_raw_points = {
        'undirected': [], 'in-in': [], 'in-out': [], 'out-in': [], 'out-out': []
    }
    
    for i in tqdm(range(num_models), desc="Analyzing Null Models"):
        file_path = os.path.join(null_dir, f"null_model_{str(i).zfill(4)}.graphml")
        if not os.path.exists(file_path):
            continue
            
        try:
            G_null = nx.read_graphml(file_path, node_type=int)
            if G_null.number_of_nodes() == 0: continue

            # Undirected case
            ann_undir = nx.average_neighbor_degree(G_null.to_undirected())
            for n, k in G_null.degree():
                if k > 0: null_ann_raw_points['undirected'].append((k, ann_undir[n]))
            
            # Directed cases
            ann_inin = nx.average_neighbor_degree(G_null, source='in', target='in')
            ann_inout = nx.average_neighbor_degree(G_null, source='in', target='out')
            ann_outin = nx.average_neighbor_degree(G_null, source='out', target='in')
            ann_outout = nx.average_neighbor_degree(G_null, source='out', target='out')

            for n, k_in in G_null.in_degree():
                if k_in > 0:
                    null_ann_raw_points['in-in'].append((k_in, ann_inin[n]))
                    null_ann_raw_points['in-out'].append((k_in, ann_inout[n]))
            for n, k_out in G_null.out_degree():
                if k_out > 0:
                    null_ann_raw_points['out-in'].append((k_out, ann_outin[n]))
                    null_ann_raw_points['out-out'].append((k_out, ann_outout[n]))

        except Exception as e:
            tqdm.write(f"Warning: Could not process {os.path.basename(file_path)}. Error: {e}")
            
    return null_ann_raw_points


def plot_ann_distribution(ax, real_data, null_data, title, x_label, y_label, y_scale='linear'):
    """
    Creates one ANN(k) vs. k plot comparing real and null model data.
    """
    # --- Plot real network data ---
    if real_data:
        k_real, ann_real = zip(*real_data)
        k_real_binned, ann_real_binned, _ = exponential_binning(real_data)
        
        ax.scatter(k_real, ann_real, color='#8F9058', s=20, alpha=0.1, label='Real Network (Raw)')
        ax.plot(k_real_binned, ann_real_binned, 's-', color='#bcbd22', markersize=8, label='Real Network (Binned)')
        
    # --- Plot null model data ---
    if null_data:
        k_null_binned, ann_null_binned, ann_null_std = exponential_binning(null_data)
        ax.errorbar(k_null_binned, ann_null_binned, yerr=ann_null_std, fmt='o',
                    color='black', markersize=5, alpha=0.8, label='Null Model (Binned)')

    # --- Formatting ---
    ax.set_xscale('log')
    ax.set_yscale(y_scale)
    ax.set_title(title, fontsize=20, fontweight='bold')
    ax.set_xlabel(x_label, fontsize=16)
    ax.set_ylabel(y_label, fontsize=16)
    ax.tick_params(labelsize=14)
    ax.legend(fontsize=12)
    ax.grid(False)


# =============================================================================
# MAIN EXECUTION
# =============================================================================

if __name__ == "__main__":
    # --- 1. Load and analyze the real network ---
    G_real = load_real_network()
    
    real_ann_data = {
        'undirected': [(k, v) for k, v in nx.average_neighbor_degree(G_real.to_undirected()).items() if k > 0],
        'in-in':      [(k, v) for k, v in nx.average_neighbor_degree(G_real, source='in', target='in').items() if k > 0],
        'in-out':     [(k, v) for k, v in nx.average_neighbor_degree(G_real, source='in', target='out').items() if k > 0],
        'out-in':     [(k, v) for k, v in nx.average_neighbor_degree(G_real, source='out', target='in').items() if k > 0],
        'out-out':    [(k, v) for k, v in nx.average_neighbor_degree(G_real, source='out', target='out').items() if k > 0]
    }
    
    # --- 2. Load and analyze null models ---
    null_ann_points = analyze_null_models_for_ann(NULL_MODELS_DIR, N_NULL_MODELS)
    
    # --- 3. Generate plots ---
    print("\nGenerating Assortativity (ANN) plots...")
    
    plot_configs = [
        {"case": "undirected", "title": "Undirected Assortativity", "xlabel": "Degree k", "ylabel": "Avg. Neighbor Degree", "yscale": "log"},
        {"case": "in-in", "title": "Directed Assortativity (In-In)", "xlabel": "In-Degree k_in", "ylabel": "Avg. Neighbor In-Degree", "yscale": "linear"},
        {"case": "in-out", "title": "Directed Assortativity (In-Out)", "xlabel": "In-Degree k_in", "ylabel": "Avg. Neighbor Out-Degree", "yscale": "linear"},
        {"case": "out-in", "title": "Directed Assortativity (Out-In)", "xlabel": "Out-Degree k_out", "ylabel": "Avg. Neighbor In-Degree", "yscale": "linear"},
        {"case": "out-out", "title": "Directed Assortativity (Out-Out)", "xlabel": "Out-Degree k_out", "ylabel": "Avg. Neighbor Out-Degree", "yscale": "log"}
    ]

    for config in plot_configs:
        fig, ax = plt.subplots(figsize=(10, 8))
        plot_ann_distribution(
            ax,
            real_ann_data[config["case"]],
            null_ann_points[config["case"]],
            config["title"],
            config["xlabel"],
            config["ylabel"],
            config["yscale"]
        )
        plt.tight_layout()
        plt.show()

    print("\nAnalysis complete.")