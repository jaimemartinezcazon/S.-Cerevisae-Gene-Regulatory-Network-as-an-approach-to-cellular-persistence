'''
Author: Jaime Martínez Cazón

Description:
This script performs a crucial data pre-processing step by converting a gene
regulatory network, represented as an adjacency matrix, into an edge list format.
The input is a tab-separated file where rows represent target genes and columns
represent regulating transcription factors. The script creates two primary
outputs:
1.  An 'edge_list.csv' file, which is the standard input for most network
    analysis tools.
2.  A 'nodes_id.csv' file, which maps each gene name to a unique integer ID.

Additionally, it performs a basic analysis to count the number of nodes that
act as both regulators and targets.
'''

import os
import pandas as pd


# =============================================================================
# SETUP: FILE PATHS
# =============================================================================
DATA_DIR = "data"
ADJACENCY_MATRIX_PATH = os.path.join(DATA_DIR, "signed_network.tsv")
OUTPUT_NODES_PATH = os.path.join(DATA_DIR, "nodes_id.csv")
OUTPUT_EDGES_PATH = os.path.join(DATA_DIR, "edge_list.csv")

# =============================================================================
# DATA CONVERSION AND ANALYSIS FUNCTIONS
# =============================================================================

def convert_adjacency_matrix_to_edgelist(df):
    """
    Converts a DataFrame representing an adjacency matrix into an edge list
    and a node ID mapping.

    Args:
        df (pd.DataFrame): The input DataFrame where rows are targets and
                           columns are sources (regulators).

    Returns:
        tuple: A tuple containing (nodes_df, edges_df).
    """
    print("Converting adjacency matrix to edge list format...")
    
    # Identify all unique nodes from both rows (targets) and columns (sources)
    source_nodes = set(df.columns)
    target_nodes = set(df.index)
    all_unique_nodes = sorted(list(source_nodes.union(target_nodes)))
    
    # Create a mapping from gene name to a unique integer ID (0-based)
    node_to_id = {node: i for i, node in enumerate(all_unique_nodes)}
    
    # Create the node ID DataFrame
    nodes_df = pd.DataFrame({
        "ID": list(node_to_id.values()),
        "Gene": list(node_to_id.keys())
    })
    
    # Create the edge list by iterating through the matrix
    # An edge exists from a column (source) to a row (target)
    edge_list = []
    for source_gene in df.columns:
        for target_gene in df.index:
            weight = df.at[target_gene, source_gene]
            if weight != 0:
                # Only include actual interactions
                edge_list.append([
                    node_to_id[source_gene], 
                    node_to_id[target_gene], 
                    int(weight)
                ])
                
    edges_df = pd.DataFrame(edge_list, columns=["Node 1", "Node 2", "Weight"])
    
    print("Conversion complete.")
    return nodes_df, edges_df


def analyze_node_roles(df):
    """
    Analyzes the roles of nodes to find those acting as both sources and targets.
    """
    print("\nAnalyzing node roles...")
    source_nodes = set(df.columns)
    target_nodes = set(df.index)
    
    common_nodes = source_nodes.intersection(target_nodes)
    num_common = len(common_nodes)
    
    print(f"Number of nodes acting as both regulators and targets: {num_common}")
    return common_nodes

# =============================================================================
# MAIN EXECUTION
# =============================================================================

if __name__ == "__main__":
    # 1. Load the adjacency matrix from the TSV file
    try:
        print(f"Loading adjacency matrix from: {ADJACENCY_MATRIX_PATH}")
        adjacency_df = pd.read_csv(ADJACENCY_MATRIX_PATH, sep="\t", index_col=0)
    except FileNotFoundError:
        print(f"Error: Input file not found at '{ADJACENCY_MATRIX_PATH}'.")
        print("Please ensure the 'signed_network.tsv' file is in the 'data' directory.")
        exit()
    except Exception as e:
        print(f"An unexpected error occurred while loading the file: {e}")
        exit()

    # 2. Analyze node roles (sources vs. targets)
    analyze_node_roles(adjacency_df)

    # 3. Convert the matrix to node and edge lists
    nodes_id_df, edge_list_df = convert_adjacency_matrix_to_edgelist(adjacency_df)
    
    # 4. Save the results to CSV files
    try:
        os.makedirs(DATA_DIR, exist_ok=True)
        nodes_id_df.to_csv(OUTPUT_NODES_PATH, index=False)
        edge_list_df.to_csv(OUTPUT_EDGES_PATH, index=False)
        
        print(f"\nSuccessfully saved node list to: '{OUTPUT_NODES_PATH}'")
        print(f"Successfully saved edge list to: '{OUTPUT_EDGES_PATH}'")
    except Exception as e:
        print(f"\nAn error occurred while saving the output files: {e}")