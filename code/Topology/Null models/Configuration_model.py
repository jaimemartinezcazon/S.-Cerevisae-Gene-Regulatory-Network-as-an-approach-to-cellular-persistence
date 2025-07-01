'''
Author: Jaime Martínez Cazón

Description:
This script generates an ensemble of surrogate networks (null models) based on 
the S. cerevisiae Gene Regulatory Network. It uses the directed configuration 
model, which preserves the in-degree and out-degree sequence of every node 
from the original network. This process randomizes the network's connections 
while maintaining its fundamental degree structure. The generated networks are 
saved in GraphML format for subsequent analysis.

Note: It is necessary to generate these null models first, before performing 
      any surrogate data analysis (the data is not uploaded because it is 
      large.)
'''

import os
import numpy as np
import pandas as pd
import networkx as nx
from tqdm import tqdm

# =============================================================================
# SETUP: CONFIGURATION PARAMETERS
# =============================================================================
# Input file path for the real network's edge list
EDGE_LIST_PATH = os.path.join("data", "giantC_edge_list.csv")

# Output directory for the generated null models
OUTPUT_DIR = "Null_Models"
NUM_NULL_MODELS = 1000  # Total number of surrogate networks to generate

# =============================================================================
# DATA LOADING FUNCTION
# =============================================================================

def load_real_network(filepath):
    """
    Loads the real network from an edge list file.

    Args:
        filepath (str): The path to the edge list CSV file.

    Returns:
        nx.DiGraph: The loaded directed graph.
    """
    if not os.path.exists(filepath):
        raise FileNotFoundError(f"The specified edge list file was not found at: {filepath}")
    
    print(f"Loading real network from: {filepath}")
    edge_list = pd.read_csv(filepath)
    
    # Create a directed graph from the edge list
    G_real = nx.from_pandas_edgelist(
        edge_list,
        source='Node 1',
        target='Node 2',
        edge_attr='Weight',
        create_using=nx.DiGraph()
    )
    print("Real network loaded successfully.")
    return G_real

# =============================================================================
# NULL MODEL GENERATION FUNCTION
# =============================================================================

def generate_null_models(G_real, num_models, output_dir):
    """
    Generates and saves a specified number of null models using the directed
    configuration model.

    Args:
        G_real (nx.DiGraph): The original network to base the models on.
        num_models (int): The number of null models to generate.
        output_dir (str): The directory where the models will be saved.
    """
    # Ensure the output directory exists
    os.makedirs(output_dir, exist_ok=True)
    print(f"Generating {num_models} null models in directory: '{output_dir}'")
    
    # Extract the degree sequences from the real network
    in_degree_sequence = [d for _, d in G_real.in_degree()]
    out_degree_sequence = [d for _, d in G_real.out_degree()]
    
    # Extract original weights to be randomly reassigned later
    original_weights = list(nx.get_edge_attributes(G_real, 'Weight').values())

    # Use tqdm for a progress bar
    for i in tqdm(range(num_models), desc="Generating Null Models"):
        # Generate a null model preserving in- and out-degree sequences.
        # This can create parallel edges and self-loops, so we use a MultiDiGraph.
        G_null_multi = nx.directed_configuration_model(
            in_degree_sequence,
            out_degree_sequence,
            create_using=nx.MultiDiGraph
        )

        # Convert to a simple DiGraph to remove parallel edges and self-loops.
        # This is a common simplification for this type of analysis.
        G_null = nx.DiGraph(G_null_multi)
        G_null.remove_edges_from(nx.selfloop_edges(G_null))

        # Randomly assign weights from the original network to the new edges.
        # We shuffle the list of original weights.
        np.random.shuffle(original_weights)

        # Assign weights to the edges in the null model.
        edges_null = list(G_null.edges())
        num_edges_to_weight = min(len(edges_null), len(original_weights))
        
        for j in range(num_edges_to_weight):
            u, v = edges_null[j]
            G_null.edges[u, v]['Weight'] = original_weights[j]

        # Save the generated null model as a GraphML file.
        # Using zfill pads the number with leading zeros (e.g., 0001, 0002).
        file_path = os.path.join(output_dir, f"null_model_{str(i).zfill(4)}.graphml")
        nx.write_graphml(G_null, file_path)

    print(f"\nSuccessfully generated and saved {num_models} null models.")

# =============================================================================
# SCRIPT EXECUTION
# =============================================================================

if __name__ == "__main__":
    # 1. Load the real network to get its properties
    try:
        real_network = load_real_network(EDGE_LIST_PATH)
        
        # 2. Generate the null models based on the real network
        generate_null_models(real_network, NUM_NULL_MODELS, OUTPUT_DIR)
        
    except FileNotFoundError as e:
        print(f"Error: {e}")
        print("Please ensure the data file is in the correct location and try again.")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")