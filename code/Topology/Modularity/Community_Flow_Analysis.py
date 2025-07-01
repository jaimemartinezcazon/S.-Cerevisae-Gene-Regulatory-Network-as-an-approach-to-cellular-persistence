'''
Author: Jaime Martínez Cazón

Description:
This script analyzes the directional nature of inter-community connections in the
S. cerevisiae Gene Regulatory Network. After detecting communities using the
undirected Louvain algorithm, it calculates the "outward flow percentage" for
each community. This metric quantifies the proportion of a community's external
connections that are outgoing, providing insight into whether a community
primarily acts as a source or a sink of regulatory information within the
network.
'''

import os
import pandas as pd
import networkx as nx
import community as community_louvain

# =============================================================================
# SETUP: FILE PATHS
# =============================================================================
DATA_DIR = "data"
EDGE_LIST_PATH = os.path.join(DATA_DIR, "giantC_edge_list.csv")

# =============================================================================
# ANALYSIS FUNCTIONS
# =============================================================================

def load_and_partition_network(edge_path):
    """
    Loads the network from an edge list, detects communities on its undirected
    version, and returns the directed graph along with the community partition.

    Args:
        edge_path (str): The path to the edge list CSV file.

    Returns:
        tuple: A tuple containing (G_directed, communities), where communities
               is a list of sets, with each set containing the nodes of a community.
    """
    print("Step 1: Loading data and detecting communities...")
    edge_list = pd.read_csv(edge_path)
    G_directed = nx.from_pandas_edgelist(
        edge_list,
        source='Node 1',
        target='Node 2',
        create_using=nx.DiGraph()
    )
    
    # Community detection is performed on the undirected version of the graph
    G_undirected = G_directed.to_undirected()
    partition_dict = community_louvain.best_partition(G_undirected, random_state=42)
    
    # Convert partition dictionary to a list of sets
    num_communities = len(set(partition_dict.values()))
    communities = [set() for _ in range(num_communities)]
    for node, cid in partition_dict.items():
        communities[cid].add(node)
        
    print(f"Found {num_communities} communities.")
    return G_directed, communities

def analyze_outward_flow(G_directed, communities):
    """
    Calculates the outward flow percentage for each detected community.

    Args:
        G_directed (nx.DiGraph): The original directed graph.
        communities (list of sets): The list of node sets for each community.

    Returns:
        pd.DataFrame: A DataFrame with the outward flow percentage for each community.
    """
    print("\nStep 2: Calculating outward flow percentage for each community...")
    
    # Using a quotient graph is an efficient way to count inter-community edges.
    # Each node in the meta-graph represents a community (a set of original nodes).
    G_meta = nx.quotient_graph(G_directed, communities, create_using=nx.MultiDiGraph)
    
    analysis_results = []
    for cid, node_set in enumerate(communities):
        # The node in the meta-graph corresponding to our community is the frozenset of its nodes.
        meta_node = frozenset(node_set)
        
        # Calculate total outgoing and incoming links from this community to others.
        # We use a MultiDiGraph to correctly sum the weights of parallel edges between communities.
        total_out_links = sum(d.get('weight', 1) for _, _, d in G_meta.out_edges(meta_node, data=True))
        total_in_links = sum(d.get('weight', 1) for _, _, d in G_meta.in_edges(meta_node, data=True))
        
        total_external_links = total_out_links + total_in_links
        
        # Calculate the outward flow percentage
        if total_external_links > 0:
            perc_outward_flow = (total_out_links / total_external_links) * 100
        else:
            # If there are no external links, the flow is undefined, but 0% is a reasonable representation.
            perc_outward_flow = 0.0
            
        analysis_results.append({
            "Community_ID": cid,
            "Total_Nodes": len(node_set),
            "Outward_Flow_Percentage": perc_outward_flow,
            "External_Outgoing_Links": total_out_links,
            "External_Incoming_Links": total_in_links,
        })
        
    return pd.DataFrame(analysis_results)

# =============================================================================
# MAIN EXECUTION
# =============================================================================

if __name__ == "__main__":
    # 1. Load the network and detect communities
    directed_graph, community_list = load_and_partition_network(EDGE_LIST_PATH)
    
    # 2. Analyze the outward flow for each community
    if directed_graph:
        results_df = analyze_outward_flow(directed_graph, community_list)
        
        # 3. Present the results in a formatted table
        print("\n--- OUTWARD FLOW PERCENTAGE BY COMMUNITY ---")
        print("(Percentage of external connections that are outgoing)")
        
        # Sort the DataFrame for better presentation
        results_df_sorted = results_df.sort_values(by="Outward_Flow_Percentage", ascending=False).reset_index(drop=True)
        
        # Format the percentage column for display
        results_df_sorted['Outward_Flow_Percentage'] = results_df_sorted['Outward_Flow_Percentage'].map('{:,.2f}%'.format)
        
        print(results_df_sorted.to_string(index=False))