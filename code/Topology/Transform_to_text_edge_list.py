'''
Author: Jaime Martínez Cazón

Description:
This script serves as a data format converter. It reads the GraphML files 
generated for the null models (which are 0-indexed and suitable for Python) 
and converts them into simple text-based edge lists. These text files follow a 
"source target weight" format with 1-based indexing, making them compatible
with the legacy Fortran code used for triangle motif analysis.
'''

import os
import networkx as nx
from tqdm import tqdm

# =============================================================================
# SETUP: CONFIGURATION PARAMETERS
# =============================================================================
# Directory where the GraphML null model files are stored
GRAPHML_DIR = "Null_Models"

# Directory where the converted .txt edge list files will be saved
EDGE_LIST_DIR = "Null_Models_EdgeLists"

# Total number of null models to process
NUM_NULL_MODELS = 1000

# =============================================================================
# CONVERSION FUNCTION
# =============================================================================

def convert_graphml_to_edgelist(graphml_dir, edgelist_dir, num_models):
    """
    Reads GraphML files, converts them to edge lists, and saves them as .txt files.

    Args:
        graphml_dir (str): The source directory containing GraphML files.
        edgelist_dir (str): The target directory for the output .txt files.
        num_models (int): The total number of models to convert.
    """
    # Ensure the output directory exists
    os.makedirs(edgelist_dir, exist_ok=True)
    
    print(f"Starting conversion of {num_models} GraphML files...")
    print(f"Source:      '{graphml_dir}'")
    print(f"Destination: '{edgelist_dir}'")
    
    # Use tqdm for a descriptive progress bar
    for i in tqdm(range(num_models), desc="Converting files"):
        # Construct the file paths for input and output
        model_number_str = str(i).zfill(4)
        graphml_path = os.path.join(graphml_dir, f"null_model_{model_number_str}.graphml")
        edgelist_path = os.path.join(edgelist_dir, f"null_model_{model_number_str}.txt")

        # Check if the source GraphML file exists
        if not os.path.exists(graphml_path):
            tqdm.write(f"Warning: File not found, skipping: {graphml_path}")
            continue

        try:
            # Read the GraphML file
            # NetworkX reads node IDs as strings from GraphML, so we'll convert them
            G = nx.read_graphml(graphml_path)

            # Open the output .txt file for writing
            with open(edgelist_path, 'w') as f:
                # Iterate through all edges and their data
                for u, v, data in G.edges(data=True):
                    # Get the weight, default to 1.0 if not present
                    weight = data.get('Weight', 1.0)
                    
                    # Convert node IDs to integers and add 1 for Fortran's 1-based indexing
                    # This is the critical step for Fortran compatibility.
                    u_fortran = int(u) + 1
                    v_fortran = int(v) + 1
                    
                    # Write the formatted line to the text file
                    f.write(f"{u_fortran} {v_fortran} {weight}\n")

        except Exception as e:
            # Report any errors during file processing
            tqdm.write(f"Error processing file {os.path.basename(graphml_path)}: {e}")

    print(f"\nConversion process completed for {num_models} models.")

# =============================================================================
# SCRIPT EXECUTION
# =============================================================================

if __name__ == "__main__":
    # Execute the conversion process
    convert_graphml_to_edgelist(GRAPHML_DIR, EDGE_LIST_DIR, NUM_NULL_MODELS)