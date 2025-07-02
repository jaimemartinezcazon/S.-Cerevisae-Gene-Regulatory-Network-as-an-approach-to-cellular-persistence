'''
Author: Jaime Martínez Cazón

Description:
This script performs a comprehensive centrality analysis on the S. cerevisiae
Gene Regulatory Network. It calculates several key centrality measures:
- Degree Centrality (in- and out-degree)
- Betweenness Centrality
- Closeness Centrality (in- and out-closeness)
- PageRank Centrality

The script then identifies the top N most central nodes for each measure and
cross-references them with functional annotations, including a list of 
candidate genes, known transcription factors (TFs), and their location within 
the network's bow-tie structure (IN, OUT, SCC, OTHERS).
'''

import os
import pandas as pd
import numpy as np
import networkx as nx
from pathlib import Path

# =============================================================================
# SETUP: FILE PATHS AND PARAMETERS
# =============================================================================
BASE_DATA_PATH = Path("data")
OUTPUT_DIR = Path("centrality_output")
N_TOP = 10  # Number of top nodes to analyze for each metric

# =============================================================================
# DATA LOADING AND PREPARATION
# =============================================================================

def load_network_and_annotations():
    """
    Loads the main network graph and all annotation files (TFs, candidates, bow-tie).
    
    Returns:
        tuple: A tuple containing the graph and various data lookups.
    """
    print("Loading all required data...")
    # Create the graph
    edge_list_path = BASE_DATA_PATH / "giantC_edge_list.csv"
    edge_list = pd.read_csv(edge_list_path)
    G = nx.from_pandas_edgelist(edge_list, 'Node 1', 'Node 2', create_using=nx.DiGraph())
    
    # Node ID to gene name mapping
    nodes_id_path = BASE_DATA_PATH / "giantC_nodes_id.csv"
    nodes_id_df = pd.read_csv(nodes_id_path)
    id_to_name = dict(zip(nodes_id_df["ID"], nodes_id_df["Gene"]))

    # Candidate genes
    try:
        candidates_path = BASE_DATA_PATH / "Candidates_list.csv" # Adjusted filename
        candidates_genes = set(pd.read_csv(candidates_path)["Gene"].dropna())
    except FileNotFoundError:
        print("Warning: Candidate list file not found. Skipping candidate analysis.")
        candidates_genes = set()

    # Transcription Factors
    try:
        tfs_path = BASE_DATA_PATH / "TFs_network.csv"
        tfs_ids = set(pd.read_csv(tfs_path)["ID"].dropna())
    except FileNotFoundError:
        print("Warning: TF list file not found. Skipping TF analysis.")
        tfs_ids = set()

    # Bow-tie sector mapping
    node_sector_map = {}
    for sector in ['IN', 'SCC', 'OUT', 'OTHERS']:
        try:
            sector_path = BASE_DATA_PATH / "Bow-Tie" / f"{sector.lower()}_sector_nodes.csv"
            sector_df = pd.read_csv(sector_path)
            for node_id in sector_df["ID"]:
                node_sector_map[node_id] = sector
        except FileNotFoundError:
            print(f"Warning: Bow-tie file for '{sector}' not found.")

    print("Data loading complete.")
    return G, nodes_id_df, id_to_name, candidates_genes, tfs_ids, node_sector_map

# =============================================================================
# CENTRALITY CALCULATION
# =============================================================================

def calculate_centrality_measures(G, nodes_df, id_to_name):
    """
    Calculates all centrality measures and returns them in a single DataFrame.
    
    Args:
        G (nx.DiGraph): The network graph.
        nodes_df (pd.DataFrame): DataFrame with node IDs and gene names.
        id_to_name (dict): Mapping from node ID to gene name.

    Returns:
        pd.DataFrame: A DataFrame containing all centrality scores for each node.
    """
    print("\nCalculating centrality measures...")
    centrality_df = nodes_df.copy()

    # All centrality calculations are unweighted
    print("  - Degree (In/Out)...")
    centrality_df["In_Degree"] = centrality_df["ID"].map(dict(G.in_degree())).fillna(0).astype(int)
    centrality_df["Out_Degree"] = centrality_df["ID"].map(dict(G.out_degree())).fillna(0).astype(int)

    print("  - Betweenness Centrality...")
    betweenness = nx.betweenness_centrality(G, weight=None)
    centrality_df["Betweenness"] = centrality_df["ID"].map(betweenness).fillna(0)

    print("  - Closeness Centrality (In/Out)...")
    closeness_out = nx.closeness_centrality(G)
    closeness_in = nx.closeness_centrality(G.reverse()) # In-closeness on reversed graph
    centrality_df["Closeness_Out"] = centrality_df["ID"].map(closeness_out).fillna(0)
    centrality_df["Closeness_In"] = centrality_df["ID"].map(closeness_in).fillna(0)

    print("  - PageRank Centrality...")
    pagerank = nx.pagerank(G, weight=None)
    centrality_df["PageRank"] = centrality_df["ID"].map(pagerank).fillna(0)
    
    # Ensure Gene_Name column exists
    centrality_df["Gene_Name"] = centrality_df["ID"].map(id_to_name)
    
    print("Centrality calculations complete.")
    return centrality_df

# =============================================================================
# TOP NODE ANALYSIS
# =============================================================================

def analyze_top_nodes(centrality_df, measure, annotations):
    """
    Analyzes and prints a report for the top N nodes for a given centrality measure.
    """
    candidates_genes, tfs_ids, node_sector_map = annotations
    
    print(f"\n--- Analysis for Top {N_TOP} nodes by {measure} ---")
    
    top_nodes = centrality_df.sort_values(by=measure, ascending=False).head(N_TOP)
    
    print(f"{'ID':<6} | {'Gene Name':<12} | {measure:<15} | {'Sector':<8} | {'Candidate':<9} | {'TF':<3}")
    print("-" * 75)
    
    sector_counts = defaultdict(int)
    
    for _, row in top_nodes.iterrows():
        node_id = row["ID"]
        gene_name = row["Gene_Name"]
        score = row[measure]
        
        sector = node_sector_map.get(node_id, "Unknown")
        is_candidate = "Yes" if gene_name in candidates_genes else "No"
        is_tf = "Yes" if node_id in tfs_ids else "No"
        
        sector_counts[sector] += 1
        
        score_str = f"{score:.6f}" if isinstance(score, float) else str(score)
        print(f"{node_id:<6} | {gene_name:<12} | {score_str:<15} | {sector:<8} | {is_candidate:<9} | {is_tf:<3}")
        
    print("-" * 75)
    print("Summary for this group:")
    for sector, count in sorted(sector_counts.items()):
        print(f"  - Sector '{sector}': {count}/{N_TOP} nodes")

# =============================================================================
# MAIN EXECUTION
# =============================================================================

if __name__ == "__main__":
    # 1. Load all data
    try:
        G_main, nodes_df_main, id_map, candidates, tfs, sector_map = load_network_and_annotations()
        annotations_tuple = (candidates, tfs, sector_map)
    except FileNotFoundError as e:
        print(f"Fatal Error: A required data file was not found. {e}")
        exit()

    # 2. Calculate all centrality measures
    centrality_results_df = calculate_centrality_measures(G_main, nodes_df_main, id_map)
    
    # 3. Save the full centrality data to a CSV
    OUTPUT_DIR.mkdir(exist_ok=True)
    output_path = OUTPUT_DIR / "full_centrality_analysis.csv"
    centrality_results_df.to_csv(output_path, index=False, float_format='%.8f')
    print(f"\nFull centrality data saved to: {output_path}")

    # 4. Analyze and report on top nodes for each measure
    centrality_measures = [
        "In_Degree", "Out_Degree", "Betweenness", 
        "Closeness_In", "Closeness_Out", "PageRank"
    ]
    
    for metric in centrality_measures:
        analyze_top_nodes(centrality_results_df, metric, annotations_tuple)

    print("\nCentrality analysis complete.")