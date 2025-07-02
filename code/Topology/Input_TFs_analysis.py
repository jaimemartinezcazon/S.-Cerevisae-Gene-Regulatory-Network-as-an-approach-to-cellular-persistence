'''
Author: Jaime Martínez Cazón

Description:
This script analyzes and compares the regulatory influence of two distinct sets
of input transcription factors (TFs), defined as having k_in=0: those associated with logarithmic growth phase
and those associated with starvation conditions. For each set of TFs, it
calculates two key metrics based on their outgoing connections in the Gene
Regulatory Network:
1.  Net Influence: The sum of the weights of all outgoing edges, indicating
    the overall activating or repressing tendency of the group.
2.  Total Out-Degree: The total number of outgoing regulatory interactions,
    representing the group's overall regulatory scope.
'''

import os
import pandas as pd

# =============================================================================
# SETUP: FILE PATHS AND GENE LISTS
# =============================================================================
DATA_DIR = "data"
EDGE_LIST_PATH = os.path.join(DATA_DIR, "giantC_edge_list.csv")
NODES_ID_PATH = os.path.join(DATA_DIR, "giantC_nodes_id.csv")

# Define the sets of Transcription Factors to be analyzed
LOG_PHASE_TFS = [
    "YDL056W", "YER111C", "YGR044C", "YIL131C", "YDR146C", "YLR131C", "YNL068C", 
    "YNL199C", "YPL075W", "YLR403W", "YPR104C", "YGL035C", "YMR172W", "YMR016C", 
    "YDR310C", "YCR084C"
]

STARVATION_TFS = [
    "YMR037C", "YGL073W", "YOL116W", "YNL027W", "YOR028C", "YDR310C", 
    "YCR084C", "YDR207C", "YFL031W"
]

# =============================================================================
# DATA LOADING AND ANALYSIS FUNCTIONS
# =============================================================================

def load_data(edge_path, nodes_path):
    """
    Loads the edge list and node ID mapping files.

    Returns:
        tuple: A tuple containing (edges_df, nodes_id_df).
    """
    print("Loading network data...")
    if not os.path.exists(edge_path) or not os.path.exists(nodes_path):
        raise FileNotFoundError("One or both required data files not found.")
    
    edges_df = pd.read_csv(edge_path)
    nodes_id_df = pd.read_csv(nodes_path, dtype={"Gene": str})
    print("Data loaded successfully.")
    return edges_df, nodes_id_df


def analyze_tf_set_influence(tf_gene_list, edges_df, nodes_id_df):
    """
    Calculates the net influence and total out-degree for a given set of TFs.

    Args:
        tf_gene_list (list): A list of gene names for the TFs.
        edges_df (pd.DataFrame): The network edge list.
        nodes_id_df (pd.DataFrame): The node ID to gene name mapping.

    Returns:
        tuple: A tuple containing (net_influence, total_out_degree).
    """
    # Step 1: Map gene names to their corresponding integer IDs
    filtered_nodes_df = nodes_id_df[nodes_id_df['Gene'].isin(tf_gene_list)]
    tf_ids = filtered_nodes_df['ID'].tolist()

    # Report any genes that were not found in the mapping
    found_genes = set(filtered_nodes_df['Gene'])
    missing_genes = [g for g in tf_gene_list if g not in found_genes]
    if missing_genes:
        print(f"    -> Warning: The following genes were not found and will be ignored: {missing_genes}")

    if not tf_ids:
        print("    -> Warning: No genes from the list were found in the network. Returning (0, 0).")
        return 0, 0

    # Step 2: Filter the edge list to find all outgoing interactions from the TF set
    outgoing_edges_df = edges_df[edges_df['Node 1'].isin(tf_ids)]
    
    # Step 3: Calculate the required metrics
    # Metric 1: Net influence (sum of weights)
    net_influence = outgoing_edges_df['Weight'].sum()
    
    # Metric 2: Total out-degree (total number of outgoing interactions)
    total_out_degree = len(outgoing_edges_df)
    
    return net_influence, total_out_degree

# =============================================================================
# MAIN EXECUTION
# =============================================================================

if __name__ == "__main__":
    try:
        # 1. Load the required data files
        edges_dataframe, nodes_dataframe = load_data(EDGE_LIST_PATH, NODES_ID_PATH)

        print("\n" + "="*60)
        print("Analyzing Outgoing Influence of TF Groups...")
        print("="*60)

        # 2. Analyze the Log Phase TF set
        print("\nProcessing Set: LogPhase_TFs")
        log_phase_influence, log_phase_degree = analyze_tf_set_influence(
            LOG_PHASE_TFS, edges_dataframe, nodes_dataframe
        )
        print(f"  -> Net Influence (Sum of Weights): {log_phase_influence}")
        print(f"  -> Total Out-Degree (Interaction Count): {log_phase_degree}")

        # 3. Analyze the Starvation TF set
        print("\nProcessing Set: Starvation_TFs")
        starvation_influence, starvation_degree = analyze_tf_set_influence(
            STARVATION_TFS, edges_dataframe, nodes_dataframe
        )
        print(f"  -> Net Influence (Sum of Weights): {starvation_influence}")
        print(f"  -> Total Out-Degree (Interaction Count): {starvation_degree}")
        
        print("\n" + "="*60)
        print("Analysis complete.")

    except FileNotFoundError as e:
        print(f"\nError: {e}")
        print("Please ensure all required data files are in the 'data' directory.")
    except Exception as e:
        print(f"\nAn unexpected error occurred: {e}")
