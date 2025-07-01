'''
Author: Jaime Martínez Cazón

Description:
This script performs a comprehensive community structure analysis of the 
S. cerevisiae Gene Regulatory Network. The workflow includes:
1.  Detecting communities using the Louvain algorithm on the undirected version
    of the network.
2.  Evaluating the statistical significance of the detected community structure
    (modularity) and individual community properties (clustering, conductance)
    by comparing them against an ensemble of null models.
3.  Locating a predefined list of candidate genes within the detected communities.
4.  Mapping the distribution of each community across the network's global
    bow-tie components (IN, OUT, SCC, OTHERS).
'''

import os
import glob
import pandas as pd
import numpy as np
import networkx as nx
import community as community_louvain
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


def load_network_and_detect_communities(edge_path):
    """Loads the network and detects communities using the Louvain algorithm."""
    print("Loading network and detecting communities...")
    edge_list = pd.read_csv(edge_path)
    G_original = nx.from_pandas_edgelist(edge_list, source='Node 1', target='Node 2', create_using=nx.DiGraph())
    G_simple = G_original.to_undirected()
    
    # Use a fixed random_state for reproducible results
    partition = community_louvain.best_partition(G_simple, random_state=42)
    num_communities = len(set(partition.values()))
    print(f"Detection complete. Found {num_communities} communities.")
    return G_original, G_simple, partition


def analyze_global_modularity(G_simple, partition, null_models_dir):
    """Analyzes the statistical significance of the global modularity."""
    print("\n--- GLOBAL MODULARITY SIGNIFICANCE ANALYSIS ---")
    real_modularity = community_louvain.modularity(partition, G_simple)
    
    null_modularity_list = []
    null_model_files = glob.glob(os.path.join(null_models_dir, "*.graphml"))
    
    for file_path in tqdm(null_model_files, desc="Analyzing Global Modularity"):
        G_null = nx.read_graphml(file_path, node_type=int).to_undirected()
        partition_null = community_louvain.best_partition(G_null, random_state=42)
        null_modularity_list.append(community_louvain.modularity(partition_null, G_null))
        
    mean_modularity = np.mean(null_modularity_list)
    std_modularity = np.std(null_modularity_list)
    p_value = calculate_empirical_p_value(real_modularity, null_modularity_list, 'greater')
    
    print(f"Real Network Modularity:         {real_modularity:.4f}")
    print(f"Null Model Modularity (Mean±SD): {mean_modularity:.4f} ± {std_modularity:.4f}")
    print(f"Empirical P-value:               {p_value:.4f}")


def analyze_community_properties(G_simple, partition, null_models_dir):
    """Analyzes properties of each community against null models."""
    print("\n--- INDIVIDUAL COMMUNITY SIGNIFICANCE ANALYSIS ---")
    nodes_by_community = {cid: {n for n, c in partition.items() if c == cid} for cid in set(partition.values())}
    
    # Aggregate metrics from all null models first
    null_metrics = {cid: {'clustering': [], 'conductance': []} for cid in nodes_by_community}
    for file_path in tqdm(glob.glob(os.path.join(null_models_dir, "*.graphml")), desc="Analyzing Community Properties"):
        G_null = nx.read_graphml(file_path, node_type=int).to_undirected()
        for cid, nodes in nodes_by_community.items():
            subgraph_null = G_null.subgraph(nodes)
            try:
                null_metrics[cid]['clustering'].append(nx.average_clustering(subgraph_null))
                null_metrics[cid]['conductance'].append(nx.conductance(G_null, nodes))
            except (nx.NetworkXError, ZeroDivisionError):
                null_metrics[cid]['clustering'].append(0.0)
                null_metrics[cid]['conductance'].append(0.0)

    # Now calculate and print results for each community
    community_stats = []
    for cid, nodes in sorted(nodes_by_community.items()):
        subgraph_real = G_simple.subgraph(nodes)
        real_clustering = nx.average_clustering(subgraph_real)
        real_conductance = nx.conductance(G_simple, nodes)
        
        p_val_clust = calculate_empirical_p_value(real_clustering, null_metrics[cid]['clustering'], 'greater')
        p_val_cond = calculate_empirical_p_value(real_conductance, null_metrics[cid]['conductance'], 'less')

        stats = {
            'Community_ID': cid, 'Num_Nodes': len(nodes), 'Density': nx.density(subgraph_real),
            'Clustering_Real': real_clustering, 'Clustering_PValue': p_val_clust,
            'Conductance_Real': real_conductance, 'Conductance_PValue': p_val_cond,
        }
        community_stats.append(stats)

    summary_df = pd.DataFrame(community_stats).set_index('Community_ID')
    print("\nSummary of Community Properties and Significance:")
    print(summary_df.to_string(float_format="%.4f"))
    summary_df.to_csv("community_significance_analysis.csv")
    print("\nDetailed summary saved to 'community_significance_analysis.csv'")


def locate_candidate_genes(partition, nodes_id_path):
    """Locates a list of candidate genes within the detected communities."""
    print("\n--- LOCATING CANDIDATE GENES ---")
    candidates_orf = ["YAL049C", "YBR208C", "YDL182W", "YFL014W", "YGR088W", "YGR180C", 
                      "YHL034C", "YJR096W", "YJR137C", "YKL001C", "YLR178C", "YML128C", 
                      "YMR105C", "YPL223C", "YPL226W"]
    
    nodes_id = pd.read_csv(nodes_id_path)
    gene_to_id_map = {gene.upper(): id_val for id_val, gene in nodes_id.set_index('ID')['Gene'].items()}

    for gene in candidates_orf:
        node_id = gene_to_id_map.get(gene.upper())
        if node_id is None:
            print(f"- {gene}: Not found in network's node list.")
        else:
            community_id = partition.get(node_id)
            if community_id is not None:
                print(f"- {gene}: Found in Community {community_id}.")
            else:
                print(f"- {gene}: Found in network (ID: {node_id}), but not in a community (isolated node).")


def map_communities_to_bowtie(partition, bowtie_dir):
    """Analyzes the distribution of communities across bow-tie components."""
    print("\n--- MAPPING COMMUNITIES TO BOW-TIE COMPONENTS ---")
    bow_tie_sets = {}
    for component in ['IN', 'SCC', 'OUT', 'OTHERS']:
        filepath = os.path.join(bowtie_dir, f"{component.lower()}_sector_nodes.csv")
        if os.path.exists(filepath):
            bow_tie_sets[component] = set(pd.read_csv(filepath)['ID'])
        else:
            print(f"Warning: Bow-tie component file not found: {filepath}")
            bow_tie_sets[component] = set()

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


# =============================================================================
# MAIN EXECUTION
# =============================================================================

if __name__ == "__main__":
    # 1. Load network and detect communities
    G_main, G_simple_main, partition_main = load_network_and_detect_communities(EDGE_LIST_PATH)
    
    # 2. Analyze global modularity significance
    analyze_global_modularity(G_simple_main, partition_main, NULL_MODELS_DIR)
    
    # 3. Analyze properties of each community
    analyze_community_properties(G_simple_main, partition_main, NULL_MODELS_DIR)
    
    # 4. Locate specific genes of interest
    locate_candidate_genes(partition_main, NODES_ID_PATH)
    
    # 5. Map communities to the bow-tie structure
    map_communities_to_bowtie(partition_main, BOWTIE_DIR)
    
    print("\n\nCommunity analysis complete.")
