'''
Author: Jaime Martínez Cazón

Description:
This script generates a 3D visualization of the S. cerevisiae Gene Regulatory
Network. It highlights Transcription Factors (TFs), sizing them based on their
out-degree to represent their influence. Nodes are distributed on the surface
of a sphere using the Fibonacci sphere algorithm for a uniform and aesthetically
pleasing layout. The visualization distinguishes between activating (gray) and
inhibiting (red) connections and uses visual effects like halos to emphasize
the TFs.
'''

import os
import pandas as pd
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# =============================================================================
# SETUP: FILE PATHS
# =============================================================================
DATA_DIR = "data"
EDGE_LIST_PATH = os.path.join(DATA_DIR, "giantC_edge_list.csv")
TFS_PATH = os.path.join(DATA_DIR, "TFs_network.csv")
DEGREES_PATH = os.path.join(DATA_DIR, "nodes_degree.csv")

# =============================================================================
# DATA LOADING AND PREPARATION
# =============================================================================

def load_data_for_visualization():
    """Loads all necessary data files."""
    print("Loading data files...")
    edges_df = pd.read_csv(EDGE_LIST_PATH)
    tfs_df = pd.read_csv(TFS_PATH)
    degrees_df = pd.read_csv(DEGREES_PATH)
    return edges_df, tfs_df, degrees_df

def fibonacci_sphere(samples):
    """Generates uniformly distributed points on a sphere."""
    points = []
    phi = np.pi * (3. - np.sqrt(5.))  # Golden angle in radians
    for i in range(samples):
        y = 1 - (i / float(samples - 1)) * 2  # y goes from 1 to -1
        radius = np.sqrt(1 - y * y)
        theta = phi * i
        x = np.cos(theta) * radius
        z = np.sin(theta) * radius
        points.append((x, y, z))
    return np.array(points)

def prepare_layout_and_attributes(G, tfs_df, degrees_df):
    """Prepares 3D positions, sizes, and colors for nodes."""
    print("Preparing layout and node attributes...")
    
    # Merge TFs with their out-degrees
    tfs_degrees = tfs_df.merge(degrees_df[['ID', 'Out']], on='ID', how='left')
    tf_nodes = tfs_degrees['ID'].tolist()
    background_nodes = [n for n in G.nodes() if n not in tf_nodes]
    
    # Generate 3D positions
    pos_3d = {}
    coords_tfs = fibonacci_sphere(len(tf_nodes))
    coords_bg = fibonacci_sphere(len(background_nodes))
    for i, node in enumerate(tf_nodes):
        pos_3d[node] = coords_tfs[i]
    for i, node in enumerate(background_nodes):
        pos_3d[node] = coords_bg[i]
        
    # Prepare TF sizes and font sizes based on out-degree
    out_degrees = tfs_degrees.set_index('ID')['Out'].to_dict()
    min_deg, max_deg = min(out_degrees.values()), max(out_degrees.values())
    
    def _map_value(val, min_val, max_val, out_min, out_max):
        if max_val == min_val: return (out_min + out_max) / 2
        return out_min + ((val - min_val) / (max_val - min_val)) * (out_max - out_min)
        
    tf_sizes = {n: _map_value(out_degrees[n], min_deg, max_deg, 100, 2000) for n in tf_nodes}
    font_sizes = {n: _map_value(out_degrees[n], min_deg, max_deg, 9, 22) for n in tf_nodes}
    
    return pos_3d, tf_nodes, background_nodes, tf_sizes, font_sizes, tfs_degrees

# =============================================================================
# VISUALIZATION FUNCTION
# =============================================================================

def create_3d_visualization(G, pos_3d, tf_nodes, bg_nodes, tf_sizes, font_sizes, tfs_info):
    """Creates and displays the 3D network visualization."""
    print("Generating 3D visualization...")
    fig = plt.figure(figsize=(16, 16))
    ax = fig.add_subplot(111, projection='3d')

    # --- Prepare colors ---
    hex_color = "#17becf"
    r, g, b = tuple(int(hex_color[i:i+2], 16)/255 for i in (1, 3, 5))
    edge_rgba = (r, g, b, 1.0)
    face_rgba = ((r + 1)/2, (g + 1)/2, (b + 1)/2, 1.0) # Lighter fill color

    # --- Draw in layers for correct depth perception ---
    
    # Layer 1: Edges
    print("  - Drawing edges...")
    for u, v, d in G.edges(data=True):
        x = [pos_3d[u][0], pos_3d[v][0]]
        y = [pos_3d[u][1], pos_3d[v][1]]
        z = [pos_3d[u][2], pos_3d[v][2]]
        color = '#d62728' if d.get('weight') == -1 else 'gray'
        alpha = 1.0 if d.get('weight') == -1 else 0.4
        ax.plot(x, y, z, color=color, alpha=alpha, linewidth=0.5)

    # Layer 2: Background nodes
    print("  - Drawing background nodes...")
    bg_coords = np.array([pos_3d[n] for n in bg_nodes])
    ax.scatter(bg_coords[:, 0], bg_coords[:, 1], bg_coords[:, 2], color='black', s=10, alpha=0.3)

    # Layer 3: TF Halos (slightly larger, edge-only circles)
    print("  - Drawing TF halos and nodes...")
    tf_coords = np.array([pos_3d[n] for n in tf_nodes])
    ax.scatter(tf_coords[:, 0], tf_coords[:, 1], tf_coords[:, 2],
               s=[tf_sizes[n] * 1.2 for n in tf_nodes],
               facecolors='none', edgecolors=[edge_rgba], linewidths=1.5, depthshade=False)

    # Layer 4: TF Nodes (main nodes with fill)
    ax.scatter(tf_coords[:, 0], tf_coords[:, 1], tf_coords[:, 2],
               s=[tf_sizes[n] for n in tf_nodes],
               facecolors=[face_rgba], edgecolors=[edge_rgba], linewidths=1.5, depthshade=False)

    # Layer 5: TF Labels (drawn last to be on top)
    print("  - Drawing TF labels...")
    id_to_genename = tfs_info.set_index('ID')['Gene'].to_dict()
    for node in tf_nodes:
        x, y, z = pos_3d[node]
        ax.text(x, y, z, id_to_genename[node],
                horizontalalignment='center', verticalalignment='center_baseline',
                fontname='Arial', fontsize=font_sizes[node], color='black', zorder=1000)

    # --- Final adjustments ---
    ax.set_axis_off()
    plt.tight_layout()
    plt.show()


# =============================================================================
# MAIN EXECUTION
# =============================================================================

if __name__ == "__main__":
    # 1. Load all necessary data
    edges, tfs, degrees = load_data_for_visualization()

    # 2. Build the graph (undirected for visualization layout)
    G_vis = nx.from_pandas_edgelist(
        edges, 'Node 1', 'Node 2', edge_attr='weight', create_using=nx.Graph()
    )

    # 3. Prepare all attributes for plotting
    (pos, tf_node_list, bg_node_list, 
     tf_size_map, font_size_map, tfs_info_df) = prepare_layout_and_attributes(G_vis, tfs, degrees)
    
    # 4. Create and display the visualization
    create_3d_visualization(G_vis, pos, tf_node_list, bg_node_list, 
                            tf_size_map, font_size_map, tfs_info_df)

    print("\nVisualization script finished.")