'''
Author: Jaime Martínez Cazón

Description:
This script creates a detailed visualization of the S. cerevisiae Gene
Regulatory Network, with nodes grouped and colored by their detected community.
The layout algorithm places communities as distinct circular clusters and uses
an iterative force-based method to prevent overlaps. The visualization is
rendered in layers to clearly distinguish intra-community and inter-community
connections.
'''

import os
import pandas as pd
import numpy as np
import networkx as nx
import community as community_louvain
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Circle
from itertools import combinations
from matplotlib.colors import LinearSegmentedColormap

# =============================================================================
# SETUP: FILE PATHS AND PARAMETERS
# =============================================================================
DATA_DIR = "data"
EDGE_LIST_PATH = os.path.join(DATA_DIR, "giantC_edge_list.csv")

# Configuration for visualization
COMMUNITY_TO_EXCLUDE = 8  # Community ID to exclude from the plot

# =============================================================================
# UTILITY AND ANALYSIS FUNCTIONS
# =============================================================================

def load_and_partition_network(edge_path):
    """Loads the network and detects communities using the Louvain algorithm."""
    print("Step 1: Loading network and detecting communities...")
    edge_list = pd.read_csv(edge_path)
    G = nx.from_pandas_edgelist(
        edge_list, source='Node 1', target='Node 2', create_using=nx.Graph()
    )
    
    partition_dict = community_louvain.best_partition(G, random_state=42)
    num_communities = len(set(partition_dict.values()))
    
    # Group nodes by community ID
    communities = {i: {n for n, cid in partition_dict.items() if cid == i} 
                   for i in range(num_communities)}
                   
    print(f"Found {num_communities} communities.")
    return G, communities, partition_dict


def resolve_overlaps(positions, radii, iterations=100, padding=0.1):
    """Iteratively adjusts circle positions to resolve overlaps."""
    print("  - Resolving community cluster overlaps...")
    pos = {k: np.array(v) for k, v in positions.items()}
    
    for _ in range(iterations):
        for c1, c2 in combinations(pos.keys(), 2):
            pos1, pos2 = pos[c1], pos[c2]
            r1, r2 = radii[c1], radii[c2]
            
            delta = pos2 - pos1
            distance = np.linalg.norm(delta)
            min_distance = r1 + r2 + padding
            
            if 0 < distance < min_distance:
                overlap = min_distance - distance
                push_vector = delta / distance
                pos[c1] -= push_vector * (overlap / 2)
                pos[c2] += push_vector * (overlap / 2)
                
    return pos


def calculate_community_layout(communities):
    """
    Calculates node positions by placing them within non-overlapping community circles.
    """
    print("Step 2: Calculating network layout...")
    nodes_per_comm = {cid: len(nodes) for cid, nodes in communities.items()}
    
    # Calculate radii for community circles based on the number of nodes
    base_radius, scale_factor = 0.05, 0.08
    comm_radii = {cid: base_radius + scale_factor * np.sqrt(count) 
                  for cid, count in nodes_per_comm.items()}
    
    # Use a meta-graph to position the community centers
    meta_graph = nx.Graph()
    meta_graph.add_nodes_from(communities.keys())
    initial_pos = nx.spring_layout(meta_graph, seed=42, k=0.8, iterations=200)
    
    # Adjust community centers to avoid overlap
    comm_pos = resolve_overlaps(initial_pos, comm_radii, padding=0.2)
    
    # Position individual nodes randomly within their community circle
    node_pos = {}
    for cid, nodes in communities.items():
        center_x, center_y = comm_pos[cid]
        radius = comm_radii[cid]
        for node in nodes:
            r = radius * np.sqrt(np.random.rand())
            theta = 2 * np.pi * np.random.rand()
            node_pos[node] = np.array([center_x + r * np.cos(theta), 
                                       center_y + r * np.sin(theta)])
            
    print("Layout calculation complete.")
    return node_pos, comm_pos, comm_radii


def create_network_visualization(G, communities, partition, node_pos, comm_pos, comm_radii):
    """
    Generates and displays the final network visualization.
    """
    print("Step 3: Generating final visualization...")
    fig, ax = plt.subplots(figsize=(24, 24))
    num_communities = len(communities)

    # --- Prepare colors and edge widths ---
    paper_colors = ['#9467bd', '#bcbd22', '#ff7f0e', '#17becf', '#d62728']
    custom_cmap = LinearSegmentedColormap.from_list("paper_palette", paper_colors, N=num_communities)
    
    G_meta = nx.quotient_graph(G, communities.values())
    inter_comm_weights = [np.log1p(d.get('weight', 1)) for _, _, d in G_meta.edges(data=True)]
    max_log_weight = max(inter_comm_weights) if inter_comm_weights else 1.0
    scaled_widths = [(16.0 * w / max_log_weight) for w in inter_comm_weights]

    # --- Draw the network in layers ---
    
    # Layer 1 & 2: Community circles and internal edges
    for cid, nodes in communities.items():
        if cid == COMMUNITY_TO_EXCLUDE: continue
        color = custom_cmap(cid / (num_communities - 1))
        
        # Background circle for the community
        circle = Circle(comm_pos[cid], comm_radii[cid], color=color, alpha=0.15, zorder=0)
        ax.add_patch(circle)
        
        # Internal edges
        nx.draw_networkx_edges(G.subgraph(nodes), node_pos, ax=ax, edge_color=color, alpha=0.5, width=0.8)

    # Layer 3: Inter-community edges
    for i, (u_nodes, v_nodes) in enumerate(G_meta.edges()):
        u_id = partition[next(iter(u_nodes))]
        v_id = partition[next(iter(v_nodes))]
        if u_id == COMMUNITY_TO_EXCLUDE or v_id == COMMUNITY_TO_EXCLUDE: continue
        
        start_pos, end_pos = comm_pos[u_id], comm_pos[v_id]
        ax.plot([start_pos[0], end_pos[0]], [start_pos[1], end_pos[1]], 
                color='gray', alpha=0.2, linewidth=scaled_widths[i], zorder=1)

    # Layer 4: Nodes
    for cid, nodes in communities.items():
        if cid == COMMUNITY_TO_EXCLUDE: continue
        color = custom_cmap(cid / (num_communities - 1))
        nx.draw_networkx_nodes(G, node_pos, nodelist=list(nodes), ax=ax, node_size=15, 
                               node_color=[color], linewidths=0)

    # Layer 5: Legend
    legend_handles = [Line2D([0], [0], marker='o', color='w', label=f"Community {i} ({len(communities[i])} nodes)",
                      markerfacecolor=custom_cmap(i / (num_communities - 1)), markersize=15)
                      for i in sorted(communities.keys()) if i != COMMUNITY_TO_EXCLUDE]
    ax.legend(handles=legend_handles, title="Communities", loc="best", fontsize=14, title_fontsize=16)

    # Final plot adjustments
    ax.set_title(f'Network Visualization by Community (Excluding Community {COMMUNITY_TO_EXCLUDE})', fontsize=28, pad=20)
    ax.set_aspect('equal', adjustable='box')
    ax.set_frame_on(False)
    ax.tick_params(left=False, bottom=False, labelleft=False, labelbottom=False)
    plt.tight_layout()
    plt.show()

# =============================================================================
# MAIN EXECUTION
# =============================================================================

if __name__ == "__main__":
    # 1. Load data and detect communities
    G_main, communities_main, partition_main = load_and_partition_network(EDGE_LIST_PATH)
    
    # 2. Calculate node and community positions
    node_positions, comm_positions, comm_radii_map = calculate_community_layout(communities_main)
    
    # 3. Generate the final plot
    create_network_visualization(G_main, communities_main, partition_main, 
                                 node_positions, comm_positions, comm_radii_map)