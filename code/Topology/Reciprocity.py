'''
Author: Jaime Martínez Cazón

Description:
This script assesses the global reciprocity of the S. cerevisiae Gene
Regulatory Network. Reciprocity measures the tendency of pairs of nodes to be
mutually connected. The script calculates the reciprocity of the real network
and compares it against an ensemble of surrogate networks (null models) to
determine if the observed reciprocity is statistically significant. A two-tailed
empirical p-value is used to evaluate whether the real network's reciprocity
is significantly different (either higher or lower) than what would be
expected by chance given the network's degree sequence.
'''

import os
import pandas as pd
import numpy as np
import networkx as nx
from tqdm import tqdm

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

def calculate_two_tailed_p_value(real_value, null_distribution):
    """
    Calculates the two-tailed empirical p-value from a null distribution.
    It assesses how many null values are as or more extreme than the real value,
    where "extremity" is measured by the absolute distance from the null mean.
    """
    n_simulations = len(null_distribution)
    if n_simulations == 0:
        return np.nan

    null_array = np.array(null_distribution)
    null_mean = np.mean(null_array)
    
    # Deviation of the real value from the null mean
    observed_deviation = abs(real_value - null_mean)
    
    # Deviations of null values from their own mean
    null_deviations = abs(null_array - null_mean)
    
    # Count how many null deviations are as or more extreme than the observed one
    count_as_extreme = np.sum(null_deviations >= observed_deviation)
    
    # The +1 correction avoids p-values of 0
    p_value = (count_as_extreme + 1) / (n_simulations + 1)
    return p_value


def load_real_network(filepath):
    """Loads the real network and extracts its giant component."""
    print(f"Loading real network from '{filepath}'...")
    if not os.path.exists(filepath):
        raise FileNotFoundError(f"Real network file not found at: {filepath}")
        
    edge_list = pd.read_csv(filepath)
    G_full = nx.from_pandas_edgelist(
        edge_list, 'Node 1', 'Node 2', create_using=nx.DiGraph()
    )
    # Ensure nodes are integers for consistency
    G_full = nx.relabel_nodes(G_full, {n: int(n) for n in G_full.nodes()})
    
    # Analysis is performed on the giant component
    giant_component_nodes = max(nx.weakly_connected_components(G_full), key=len)
    G_real = G_full.subgraph(giant_component_nodes).copy()
    
    print(f"Real network (giant component) loaded: {G_real.number_of_nodes()} nodes.")
    return G_real


def analyze_null_models_reciprocity(null_dir, num_models):
    """
    Loads and analyzes the ensemble of null models to calculate reciprocity.

    Returns:
        list: A list of reciprocity values, one for each successfully loaded model.
    """
    print(f"\nAnalyzing up to {num_models} null models for reciprocity...")
    null_reciprocity_list = []
    
    for i in tqdm(range(num_models), desc="Analyzing Null Models"):
        file_path = os.path.join(null_dir, f"null_model_{str(i).zfill(4)}.graphml")
        if not os.path.exists(file_path):
            continue
            
        try:
            # It is crucial to load the models as directed graphs for reciprocity
            G_null = nx.read_graphml(file_path, node_type=int)
            if G_null.number_of_nodes() > 0:
                null_reciprocity_list.append(nx.reciprocity(G_null))
        except Exception as e:
            tqdm.write(f"Warning: Could not process {os.path.basename(file_path)}. Error: {e}")
            
    return null_reciprocity_list

# =============================================================================
# MAIN EXECUTION
# =============================================================================

if __name__ == "__main__":
    try:
        # 1. Load the real network and calculate its reciprocity
        G_real = load_real_network(EDGE_LIST_PATH)
        real_reciprocity = nx.reciprocity(G_real)
        
        # 2. Load null models and calculate their reciprocity values
        null_reciprocity_values = analyze_null_models_reciprocity(NULL_MODELS_DIR, N_NULL_MODELS)
        
        # 3. Perform statistical significance analysis
        print("\n" + "="*70)
        print("--- GLOBAL RECIPROCITY SIGNIFICANCE ANALYSIS ---")
        
        if null_reciprocity_values:
            mean_null_reciprocity = np.mean(null_reciprocity_values)
            std_null_reciprocity = np.std(null_reciprocity_values)
            
            p_value = calculate_two_tailed_p_value(real_reciprocity, null_reciprocity_values)
            
            print(f"Real Network Reciprocity:               {real_reciprocity:.6f}")
            print(f"Null Model Reciprocity (Mean ± SD):     {mean_null_reciprocity:.6f} ± {std_null_reciprocity:.6f}")
            print(f"Two-Tailed Empirical P-value:           {p_value:.6f}")
            
            # Interpretation of the result
            alpha = 0.05
            print(f"\nSignificance Level (alpha): {alpha}")
            if p_value < alpha:
                print(f"Conclusion: The observed reciprocity ({real_reciprocity:.6f}) is statistically significant (p < {alpha}).")
                print("This indicates the network's tendency for mutual connections is not random.")
            else:
                print(f"Conclusion: The observed reciprocity is NOT statistically significant (p >= {alpha}).")
                print("This suggests the level of reciprocity could be explained by the degree sequence alone.")
        else:
            print("No null model data available to perform significance analysis.")
            
        print("="*70)

    except FileNotFoundError as e:
        print(f"\nFatal Error: {e}")
        print("Please ensure all required data files are in the correct 'data' directory.")
    except Exception as e:
        print(f"\nAn unexpected error occurred: {e}")