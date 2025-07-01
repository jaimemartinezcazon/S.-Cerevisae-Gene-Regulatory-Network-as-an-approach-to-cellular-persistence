'''
Author: Jaime Martínez Cazón

Description:
This script performs a community structure analysis on the S. cerevisiae Gene
Regulatory Network, specifically using an algorithm that accounts for edge
directionality. The workflow includes:
1.  Detecting communities using the directed Louvain algorithm available in
    NetworkX.
2.  Evaluating the statistical significance of the directed modularity and
    other community properties against a null model ensemble.
3.  Locating candidate genes within the detected directed communities.
4.  Mapping the community distribution across the network's bow-tie structure.
'''

import os
import glob
import pandas as pd
import numpy as np
import networkx as nx
from tqdm import tqdm

# =============================================================================
# SETUP: FILE PATHS AND PARAMETERS
# =============================================================================
DATA_DIR = "data"
EDGE_LIST_PATH = os.path.join(DATA_DIR, "giantC_edge_list.csv")
NODES_ID_PATH = os.path.join(DATA_DIR, "giantC_nodes_id.csv")
NULL_MODELS_DIR = os.path.join(DATA_DIR, "Null Models")
BOWTIE_DIR = os.path.join(DATA_DIR, "Bow-Tie")

# =============================================================================
# UTILITY AND ANALYSIS FUNCTIONS
# =============================================================================

def calculate_empirical_p_value(real_value, null_distribution, direction='greater'):
    """Calculates the empirical p-value from a null distribution."""
    n_simulations = len(null_distribution)
    if n_simulations == 0: return np.nan
    
    null_array = np.array(null_distribution)
    count = np.sum(null_array >= real_value) if direction == 'greater' else np.sum(null_array <= real_value)
    return (count + 1) / (n_simulations + 1)

def load_network_and_detect_directed_communities(edge_path):
    """Loads the directed network and detects communities using a directed algorithm."""
    print("Loading directed network and detecting communities...")
    edge_list = pd.read_csv(edge_path)
    G = nx.from_pandas_edgelist(edge_list, source='Node 1', target='Node 2', edge_attr='Weight', create_using=nx.DiGraph())
    
    # Use the Louvain algorithm for directed graphs from NetworkX
    communities_list = nx.community.louvain_communities(G, seed=42)
    
    # Convert the list of sets to a partition dictionary for compatibility
    partition = {node: i for i, comm in enumerate(communities_list) for node in comm}
    
    print(f"Detection complete. Found {len(communities_list)} communities.")
    return G, communities_list, partition

def analyze_directed_modularity_significance(G, communities, null_models_dir):
    """Analyzes the statistical significance of directed modularity."""
    print("\n--- DIRECTED MODULARITY SIGNIFICANCE ANALYSIS ---")
    real_modularity = nx.community.modularity(G, communities)
    
    null_modularity_list = []
    for file_path in tqdm(glob.glob(os.path.join(null_models_dir, "*.graphml")), desc="Analyzing Directed Modularity"):
        G_null = nx.read_graphml(file_path, node_type=int)
        communities_null = nx.community.louvain_communities(G_null, seed=42)
        null_modularity_list.append(nx.community.modularity(G_null, communities_null))
        
    mean_modularity = np.mean(null_modularity_list)
    std_modularity = np.std(null_modularity_list)
    p_value = calculate_empirical_p_value(real_modularity, null_modularity_list, 'greater')
    
    print(f"Real Network Directed Modularity:  {real_modularity:.4f}")
    print(f"Null Model Modularity (Mean±SD):  {mean_modularity:.4f} ± {std_modularity:.4f}")
    print(f"Empirical P-value:                 {p_value:.4f}")

def analyze_directed_community_properties(G, partition, null_models_dir):
    """Analyzes directed properties of each community against null models."""
    print("\n--- INDIVIDUAL DIRECTED COMMUNITY ANALYSIS ---")
    nodes_by_community = {cid: {n for n, c in partition.items() if c == cid} for cid in set(partition.values())}
    
    null_metrics = {cid: {'clustering': [], 'conductance': []} for cid in nodes_by_community}
    for file_path in tqdm(glob.glob(os.path.join(null_models_dir, "*.graphml")), desc="Analyzing Community Properties"):
        G_null = nx.read_graphml(file_path, node_type=int)
        for cid, nodes in nodes_by_community.items():
            subgraph_null = G_null.subgraph(nodes)
            try:
                null_metrics[cid]['clustering'].append(nx.average_clustering(subgraph_null))
                null_metrics[cid]['conductance'].append(nx.conductance(G_null, nodes))
            except (nx.NetworkXError, ZeroDivisionError):
                null_metrics[cid]['clustering'].append(0.0)
                null_metrics[cid]['conductance'].append(0.0)

    community_stats = []
    for cid, nodes in sorted(nodes_by_community.items()):
        subgraph_real = G.subgraph(nodes)
        real_clustering = nx.average_clustering(subgraph_real)
        real_conductance = nx.conductance(G, nodes)
        
        p_val_clust = calculate_empirical_p_value(real_clustering, null_metrics[cid]['clustering'], 'greater')
        p_val_cond = calculate_empirical_p_value(real_conductance, null_metrics[cid]['conductance'], 'less')
        
        stats = {'Community_ID': cid, 'Num_Nodes': len(nodes), 'Density': nx.density(subgraph_real),
                 'Clustering_Real': real_clustering, 'Clustering_PValue': p_val_clust,
                 'Conductance_Real': real_conductance, 'Conductance_PValue': p_val_cond}
        community_stats.append(stats)
        
    summary_df = pd.DataFrame(community_stats).set_index('Community_ID')
    print("\nSummary of Directed Community Properties and Significance:")
    print(summary_df.to_string(float_format="%.4f"))
    summary_df.to_csv("community_significance_directed.csv")
    print("\nDetailed summary saved to 'community_significance_directed.csv'")

# Functions 'locate_candidate_genes' and 'map_communities_to_bowtie' are reused from the undirected script,
# as they only depend on the final partition, not the graph type itself.

def locate_candidate_genes(partition, nodes_id_path):
    """Locates a list of candidate genes within the detected communities."""
    # (This function is identical to the one in the undirected analysis script)
    print("\n--- LOCATING CANDIDATE GENES ---")
    candidates_orf = ["YAL049C", "YBR208C", "YDL182W", "YFL014W", "YGR088W", "YGR180C", 
                      "YHL034C", "YJR096W", "YJR137C", "YKL001C", "YLR178C", "YML128C", 
                      "YMR105C", "YPL223C", "YPL226W"]
    
    nodes_id = pd.read_csv(nodes_id_path)
    gene_to_id_map = {gene.upper(): id_val for id_val, gene in nodes_id.set_index('ID')['Gene'].items()}

    for gene in candidates_orf:
        node_id = gene_to_id_map.get(gene.upper())
        if node_id:
            community_id = partition.get(node_id, 'N/A (Isolated)')
            print(f"- {gene}: Found in Community {community_id}.")
        else:
            print(f"- {gene}: Not found in node list.")

def map_communities_to_bowtie(partition, bowtie_dir):
    """Analyzes the distribution of communities across bow-tie components."""
    # (This function is identical to the one in the undirected analysis script)
    print("\n--- MAPPING COMMUNITIES TO BOW-TIE COMPONENTS ---")
    bow_tie_sets = {comp: set(pd.read_csv(os.path.join(bowtie_dir, f"{comp.lower()}_sector_nodes.csv"))['ID'])
                    for comp in ['IN', 'SCC', 'OUT', 'OTHERS']}
    
    nodes_by_community = {cid: {n for n, c in partition.items() if c == cid} for cid in set(partition.values())}
    
    analysis_results = []
    for cid, nodes in sorted(nodes_by_community.items()):
        total_nodes = len(nodes)
        if total_nodes == 0: continue
        
        counts = {comp: len(nodes.intersection(s)) for comp, s in bow_tie_sets.items()}
        result_row = {'Community_ID': cid, 'Total_Nodes': total_nodes}
        result_row.update({f'% {comp}': (count / total_nodes) * 100 for comp, count in counts.items()})
        analysis_results.append(result_row)
        
    results_df = pd.DataFrame(analysis_results).set_index('Community_ID')
    print("\nDistribution of Communities across Bow-Tie Components (%):")
    print(results_df.to_string(float_format="%.2f"))
    results_df.to_csv("community_bowtie_distribution_directed.csv")
    print("\nBow-tie distribution report saved to 'community_bowtie_distribution_directed.csv'")

# =============================================================================
# MAIN EXECUTION
# =============================================================================

if __name__ == "__main__":
    # 1. Load directed network and detect directed communities
    G_main, communities_main, partition_main = load_network_and_detect_directed_communities(EDGE_LIST_PATH)
    
    # 2. Analyze directed modularity significance
    analyze_directed_modularity_significance(G_main, communities_main, NULL_MODELS_DIR)
    
    # 3. Analyze properties of each directed community
    analyze_directed_community_properties(G_main, partition_main, NULL_MODELS_DIR)
    
    # 4. Locate specific genes of interest
    locate_candidate_genes(partition_main, NODES_ID_PATH)
    
    # 5. Map communities to the bow-tie structure
    map_communities_to_bowtie(partition_main, BOWTIE_DIR)
    
    print("\n\nDirected community analysis complete.")