'''
Author: Jaime Martínez Cazón

Description:
This script performs a comprehensive topological analysis of the S. cerevisiae 
Gene Regulatory Network. It covers several key aspects of the network's 
structure:
1.  Global degree distributions (in-degree and out-degree) and power-law fits.
2.  Analysis of the Strongly Connected Component (SCC), including the degree
    properties of its internal and external connections.
3.  Degree distribution analysis for each of the main bow-tie components 
    (IN, OUT, SCC, OTHERS).
4.  A specific analysis of the inhibitory subnetwork.
'''

import os
import pandas as pd
import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
import powerlaw

# =============================================================================
# SETUP: FILE PATHS
# =============================================================================
DATA_DIR = "data"
EDGE_LIST_PATH = os.path.join(DATA_DIR, "giantC_edge_list.csv")
NODES_ID_PATH = os.path.join(DATA_DIR, "giantC_nodes_id.csv")
BOWTIE_COMPONENTS_DIR = os.path.join(DATA_DIR, "Bow-Tie")

# =============================================================================
# DATA LOADING
# =============================================================================

def load_network_data():
    """
    Loads the network edge list and node ID mappings.

    Returns:
        tuple: A tuple containing the NetworkX DiGraph and a dictionary
               mapping node IDs to gene names.
    """
    print("Loading network data...")
    edge_list = pd.read_csv(EDGE_LIST_PATH)
    nodes_id = pd.read_csv(NODES_ID_PATH)
    
    id_to_name = dict(zip(nodes_id["ID"], nodes_id["Gene"]))
    
    G = nx.from_pandas_edgelist(
        edge_list,
        source='Node 1',
        target='Node 2',
        edge_attr='Weight',
        create_using=nx.DiGraph()
    )
    print("Network data loaded successfully.")
    return G, id_to_name

# =============================================================================
# ANALYSIS FUNCTIONS
# =============================================================================

def analyze_global_properties(G):
    """
    Calculates and prints basic global properties of the network.
    """
    print("\n--- GLOBAL NETWORK PROPERTIES ---")
    avg_in_degree = np.mean([d for _, d in G.in_degree()])
    avg_out_degree = np.mean([d for _, d in G.out_degree()])
    print(f"Average in-degree: {avg_in_degree:.2f}")
    print(f"Average out-degree: {avg_out_degree:.2f}")
    
    num_edges_neg1 = sum(1 for _, _, d in G.edges(data=True) if d.get("Weight") == -1)
    percent_edges_neg1 = (num_edges_neg1 / G.number_of_edges()) * 100
    print(f"Percentage of inhibitory edges (weight -1): {percent_edges_neg1:.2f}%")
    print("-" * 35)

def plot_global_degree_distributions(G):
    """
    Plots the global in-degree and out-degree distributions (PDF, CCDF) and fits.
    """
    print("\nPlotting global degree distributions...")
    degree_in = dict(G.in_degree())
    degree_out = dict(G.out_degree())
    
    # --- Power-law fit ---
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    
    def _plot_fit(ax, degree_dict, degree_type):
        data = np.array([d for d in degree_dict.values() if d > 0])
        fit = powerlaw.Fit(data, verbose=False)
        fit.plot_pdf(ax=ax, marker='o', linestyle='None', label="Data")
        fit.power_law.plot_pdf(ax=ax, linestyle='-', color='r', label=f'Fit (γ={fit.alpha:.2f})')
        ax.set_title(f"Power-Law Fit ({degree_type}-Degree)", fontsize=20, fontweight='bold')
        ax.set_xlabel("Degree (k)", fontsize=16)
        ax.set_ylabel("P(k)", fontsize=16)
        ax.legend(fontsize=14)

    _plot_fit(axes[0], degree_in, "In")
    _plot_fit(axes[1], degree_out, "Out")
    plt.tight_layout()
    plt.show()

def analyze_scc_properties(G, id_to_name):
    """
    Analyzes the properties of the largest Strongly Connected Component (SCC).
    """
    print("\nAnalyzing Strongly Connected Component (SCC)...")
    largest_scc_nodes = max(nx.strongly_connected_components(G), key=len)
    G_scc = G.subgraph(largest_scc_nodes).copy()
    
    # Get degrees and sort nodes by global out-degree
    nodes_sorted = sorted(G_scc.nodes(), key=lambda n: G.out_degree(n), reverse=True)
    gene_names_sorted = [id_to_name.get(n, str(n)) for n in nodes_sorted]
    
    # Prepare data for plotting
    data = {
        'Global In': [G.in_degree(n) for n in nodes_sorted],
        'Global Out': [G.out_degree(n) for n in nodes_sorted],
        'Internal In': [G_scc.in_degree(n) for n in nodes_sorted],
        'Internal Out': [G_scc.out_degree(n) for n in nodes_sorted]
    }
    
    # Plotting function
    def _plot_scc_bars(ax, title, in_data, out_data):
        x = np.arange(len(nodes_sorted))
        width = 0.4
        ax.bar(x - width/2, in_data, width, label='In-Degree')
        ax.bar(x + width/2, out_data, width, label='Out-Degree')
        ax.set_title(title, fontsize=20, fontweight='bold')
        ax.set_ylabel("Degree", fontsize=16)
        ax.set_xticks(x)
        ax.set_xticklabels(gene_names_sorted, rotation=90, fontsize=10)
        ax.legend(fontsize=14)

    fig, axes = plt.subplots(1, 2, figsize=(18, 7), sharey=True)
    _plot_scc_bars(axes[0], "SCC Global Connections", data['Global In'], data['Global Out'])
    _plot_scc_bars(axes[1], "SCC Internal Connections", data['Internal In'], data['Internal Out'])
    plt.tight_layout()
    plt.show()

def analyze_bowtie_components(G):
    """
    Analyzes the degree distributions for each bow-tie component.
    """
    print("\nAnalyzing degree distributions per bow-tie component...")
    
    def _plot_degree_distribution_fit(ax, data, title):
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_title(title, fontsize=18, fontweight='bold')
        ax.set_xlabel("Degree (k)", fontsize=14)
        ax.set_ylabel("P(k)", fontsize=14)
        if len(data) > 1:
            fit = powerlaw.Fit(data, verbose=False)
            fit.plot_pdf(ax=ax, marker='o', linestyle='None', label='Data')
            fit.power_law.plot_pdf(ax=ax, linestyle='-', color='r', label=f'Fit (γ={fit.alpha:.2f})')
            ax.legend()
        else:
            ax.text(0.5, 0.5, "Insufficient data", ha='center', va='center', transform=ax.transAxes)

    components = {'IN': 'in_sector_nodes.csv', 'SCC': 'scc_nodes.csv',
                  'OUT': 'out_sector_nodes.csv', 'OTHERS': 'others_nodes.csv'}

    for comp_name, filename in components.items():
        filepath = os.path.join(BOWTIE_COMPONENTS_DIR, filename)
        if not os.path.exists(filepath):
            print(f"Warning: Component file not found, skipping: {filepath}")
            continue
            
        ids = pd.read_csv(filepath)['ID'].tolist()
        deg_in = np.array([G.in_degree(n) for n in ids if G.in_degree(n) > 0])
        deg_out = np.array([G.out_degree(n) for n in ids if G.out_degree(n) > 0])
        
        fig, axes = plt.subplots(1, 2, figsize=(14, 6))
        fig.suptitle(f"Degree Distribution for {comp_name} Component", fontsize=22, fontweight='bold')
        _plot_degree_distribution_fit(axes[0], deg_in, "In-Degree Distribution")
        _plot_degree_distribution_fit(axes[1], deg_out, "Out-Degree Distribution")
        plt.tight_layout(rect=[0, 0.03, 1, 0.95])
        plt.show()

# =============================================================================
# MAIN EXECUTION
# =============================================================================

if __name__ == "__main__":
    # Load the network data
    G_main, id_to_name_map = load_network_data()
    
    # Run the various analyses
    analyze_global_properties(G_main)
    plot_global_degree_distributions(G_main)
    analyze_scc_properties(G_main, id_to_name_map)
    analyze_bowtie_components(G_main)
    
    print("\nTopological analysis complete.")