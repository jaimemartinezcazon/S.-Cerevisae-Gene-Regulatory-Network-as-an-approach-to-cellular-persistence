'''
Author: Jaime Martínez Cazón

Description:
This script performs a bow-tie decomposition analysis on the S. cerevisiae 
Gene Regulatory Network. It calculates the size of the main bow-tie components 
(IN, OUT, SCC, and OTHERS) for the real network and compares these results 
against an ensemble of surrogate networks (null models) generated using a 
configuration model. The script then generates a summary table and a bar plot 
to visualize the comparison.
'''

import os
import json
import time
import numpy as np
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from collections import deque

# =============================================================================
# SETUP: DIRECTORIES AND FILE PATHS
# =============================================================================
# The code assumes a project structure where data is stored in a 'data' folder.
# All paths are relative to the project's root directory.

# Input data paths
DATA_DIR = "data"
EDGE_LIST_PATH = os.path.join(DATA_DIR, "giantC_edge_list.csv")
NODES_ID_PATH = os.path.join(DATA_DIR, "giantC_nodes_id.csv")
NULL_MODELS_DIR = os.path.join(DATA_DIR, "Null Models")

# Output file paths
RESULTS_FILE = "bow_tie_comparison_results.json"
PLOT_DATA_FILE = "bow_tie_plot_percentages.json"

# =============================================================================
# UTILITY FUNCTIONS
# =============================================================================

def multi_source_shortest_path_length(G, sources):
    """
    Computes the shortest path length from a set of source nodes to all other
    reachable nodes in an unweighted graph using Breadth-First Search (BFS).

    Args:
        G (nx.DiGraph): The graph to search.
        sources (iterable): A collection of source nodes.

    Returns:
        dict: A dictionary mapping reachable nodes to their shortest distance
              from the nearest source node.
    """
    dist = {}
    queue = deque()
    for s in sources:
        if s in G:
            dist[s] = 0
            queue.append(s)

    while queue:
        u = queue.popleft()
        for v in G.neighbors(u):
            if v not in dist:
                dist[v] = dist[u] + 1
                queue.append(v)
    return dist


def calculate_bow_tie_metrics(G):
    """
    Calculates key metrics for a bow-tie decomposition of a directed graph.
    The components are IN, OUT, SCC (Strongly Connected Component), and OTHERS.

    Args:
        G (nx.DiGraph): The directed graph to analyze.

    Returns:
        dict: A dictionary containing the calculated bow-tie metrics. Returns
              None if a meaningful analysis cannot be performed (e.g., no SCC).
    """
    N = G.number_of_nodes()
    if N == 0:
        return None

    # Identify the largest strongly connected component (SCC)
    sccs = sorted(nx.strongly_connected_components(G), key=len, reverse=True)
    if not sccs or len(sccs[0]) < 2:
        print("Warning: Graph has no significant SCC (size < 2).")
        return None

    scc_nodes_set = set(sccs[0])

    # Compute IN and OUT components based on reachability from the SCC
    # OUT: Nodes reachable from the SCC (descendants)
    reachable_from_scc = set()
    for node in scc_nodes_set:
        reachable_from_scc.update(nx.descendants(G, node))
    out_sector_set = reachable_from_scc - scc_nodes_set

    # IN: Nodes that can reach the SCC (ancestors)
    # This is efficiently computed by finding descendants in the reversed graph
    G_rev = G.reverse(copy=False)
    reachable_to_scc = set()
    for node in scc_nodes_set:
        reachable_to_scc.update(nx.descendants(G_rev, node))
    in_sector_set = reachable_to_scc - scc_nodes_set

    # OTHERS: All nodes not in IN, SCC, or OUT
    all_nodes_set = set(G.nodes())
    others_set = all_nodes_set - scc_nodes_set - in_sector_set - out_sector_set

    # --- Calculate component sizes and percentages ---
    in_size = len(in_sector_set)
    scc_size = len(scc_nodes_set)
    out_size = len(out_sector_set)
    others_size = len(others_set)
    
    # --- Calculate interface densities (edges per node) ---
    # IN -> SCC interface density
    in_to_scc_edges = sum(1 for u, v in G.edges() if u in in_sector_set and v in scc_nodes_set)
    ratio_in_scc = in_to_scc_edges / in_size if in_size > 0 else 0.0

    # SCC -> OUT interface density
    scc_to_out_edges = sum(1 for u, v in G.edges() if u in scc_nodes_set and v in out_sector_set)
    ratio_scc_out = scc_to_out_edges / out_size if out_size > 0 else 0.0

    # --- Calculate average shortest path lengths between components ---
    # IN -> SCC average distance
    dist_to_scc = multi_source_shortest_path_length(G_rev, scc_nodes_set)
    distances_in_to_scc = [dist_to_scc[node] for node in in_sector_set if node in dist_to_scc]
    avg_dist_in_scc = np.mean(distances_in_to_scc) if distances_in_to_scc else float('nan')

    # SCC -> OUT average distance
    dist_from_scc = multi_source_shortest_path_length(G, scc_nodes_set)
    distances_scc_to_out = [dist_from_scc[node] for node in out_sector_set if node in dist_from_scc]
    avg_dist_scc_out = np.mean(distances_scc_to_out) if distances_scc_to_out else float('nan')

    # --- Compile all metrics into a dictionary ---
    metrics = {
        'N_nodes': N,
        'L_edges': G.number_of_edges(),
        'in_pct': (in_size / N) * 100,
        'scc_pct': (scc_size / N) * 100,
        'out_pct': (out_size / N) * 100,
        'others_pct': (others_size / N) * 100,
        'ratio_in_scc': ratio_in_scc,
        'ratio_scc_out': ratio_scc_out,
        'average_dist_in_scc': avg_dist_in_scc,
        'average_dist_scc_out': avg_dist_scc_out
    }
    return metrics

# =============================================================================
# MAIN ANALYSIS SCRIPT
# =============================================================================

def main():
    """
    Main function to execute the full analysis pipeline.
    """
    start_time = time.time()
    
    # --- 1. Load and build the real network ---
    print("Loading real network...")
    edge_list = pd.read_csv(EDGE_LIST_PATH)
    nodes_id_df = pd.read_csv(NODES_ID_PATH)
    
    G_real = nx.DiGraph()
    G_real.add_nodes_from(nodes_id_df['ID'].astype(int))
    for _, row in edge_list.iterrows():
        G_real.add_edge(int(row['Node 1']), int(row['Node 2']), weight=int(row['Weight']))
    
    print(f"Real network loaded. Nodes: {G_real.number_of_nodes()}, Edges: {G_real.number_of_edges()}")

    # --- 2. Check for pre-computed results or compute new ones ---
    if os.path.exists(RESULTS_FILE):
        print(f"Loading previously computed results from {RESULTS_FILE}")
        with open(RESULTS_FILE, 'r') as f:
            results_data = json.load(f)
        real_metrics = results_data['real_metrics']
        null_mean_metrics = results_data['null_mean_metrics']
        null_std_metrics = results_data['null_std_metrics']
    else:
        print("\nComputing metrics for the real network...")
        real_metrics = calculate_bow_tie_metrics(G_real)

        print("Computing metrics for null models...")
        null_metrics_list = []
        null_files = [f for f in os.listdir(NULL_MODELS_DIR) if f.endswith(".graphml")]
        
        for i, file_name in enumerate(null_files):
            if (i + 1) % 100 == 0:
                print(f"Processing null model {i + 1}/{len(null_files)}")
            
            file_path = os.path.join(NULL_MODELS_DIR, file_name)
            G_null = nx.read_graphml(file_path)
            G_null = nx.relabel_nodes(G_null, {n: int(n) for n in G_null.nodes()})
            
            metrics = calculate_bow_tie_metrics(G_null)
            if metrics:
                null_metrics_list.append(metrics)
        
        # Calculate mean and standard deviation for the null model ensemble
        if not null_metrics_list:
            raise ValueError("No valid null model metrics were collected.")
        
        df_null = pd.DataFrame(null_metrics_list)
        null_mean_metrics = df_null.mean().to_dict()
        null_std_metrics = df_null.std().to_dict()

        # Save computed results to a JSON file for future use
        print(f"\nSaving results to {RESULTS_FILE}...")
        results_data = {
            'real_metrics': real_metrics,
            'null_mean_metrics': null_mean_metrics,
            'null_std_metrics': null_std_metrics
        }
        with open(RESULTS_FILE, 'w') as f:
            json.dump(results_data, f, indent=4)

    # --- 3. Print a summary comparison table ---
    print("\n--- Bow-Tie Metrics Comparison ---")
    metrics_to_print = {
        'in_pct': '% IN Nodes',
        'scc_pct': '% SCC Nodes',
        'out_pct': '% OUT Nodes',
        'others_pct': '% OTHERS Nodes',
        'ratio_in_scc': 'Ratio Edges IN->SCC / IN Node',
        'ratio_scc_out': 'Ratio Edges SCC->OUT / OUT Node',
        'average_dist_in_scc': 'Avg. Dist. IN -> SCC',
        'average_dist_scc_out': 'Avg. Dist. SCC -> OUT',
    }

    print(f"{'Metric':<35} | {'Real Network':<15} | {'Null Models (Mean ± Std)'}")
    print("-" * 85)
    for key, label in metrics_to_print.items():
        real_val = real_metrics.get(key, float('nan'))
        null_mean = null_mean_metrics.get(key, float('nan'))
        null_std = null_std_metrics.get(key, float('nan'))
        print(f"{label:<35} | {real_val:<15.4f} | {null_mean:.4f} ± {null_std:.4f}")
    print("-" * 85)
    
    # --- 4. Prepare and save data required for plotting ---
    plot_sectors = ['in_pct', 'scc_pct', 'out_pct', 'others_pct']
    plot_labels = ['IN', 'SCC', 'OUT', 'Others']

    plot_data = {
        'labels': plot_labels,
        'real_percentages': [real_metrics.get(sec, 0.0) for sec in plot_sectors],
        'null_mean_percentages': [null_mean_metrics.get(sec, 0.0) for sec in plot_sectors],
        'null_std_percentages': [null_std_metrics.get(sec, 0.0) for sec in plot_sectors]
    }
    with open(PLOT_DATA_FILE, 'w') as f:
        json.dump(plot_data, f, indent=4)
    print(f"Plotting data saved to {PLOT_DATA_FILE}")

    print(f"\nAnalysis complete in {time.time() - start_time:.2f} seconds.")
    return plot_data

# =============================================================================
# PLOTTING
# =============================================================================

def create_plot(plot_data):
    """
    Generates and displays the bar plot comparing bow-tie sector percentages.
    
    Args:
        plot_data (dict): A dictionary containing the data needed for the plot.
    """
    if not plot_data:
        print("No data available for plotting.")
        return

    labels = plot_data['labels']
    real_percentages = plot_data['real_percentages']
    null_mean_percentages = plot_data['null_mean_percentages']
    null_std_percentages = plot_data['null_std_percentages']

    x = np.arange(len(labels))
    width = 0.35

    fig, ax = plt.subplots(figsize=(12, 8))

    # Bars for the real network
    ax.bar(x - width/2, real_percentages, width, label='Real Network', color='#bcbd22')

    # Bars for the null models (surrogate data)
    ax.bar(x + width/2, null_mean_percentages, width, 
           yerr=null_std_percentages, capsize=5,
           label='Surrogate Data (Mean ± Std)', color='#cccccc',
           hatch='//', edgecolor='black', alpha=0.8)

    # Plot formatting
    ax.set_ylabel('% of Nodes', fontsize=20)
    ax.set_title('Bow-Tie Sector Node Percentage Comparison', fontsize=24, fontweight='bold')
    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=18)
    ax.tick_params(axis='y', labelsize=16)
    ax.legend(fontsize=16)
    fig.tight_layout()
    plt.show()

# =============================================================================
# SCRIPT EXECUTION
# =============================================================================

if __name__ == "__main__":
    # Execute the main analysis pipeline
    plotting_data = main()
    
    # Generate the plot using the results
    create_plot(plotting_data)
