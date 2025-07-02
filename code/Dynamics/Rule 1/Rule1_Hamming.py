'''
Author: Jaime Martínez Cazón

Description:
This script performs a large-scale dynamical simulation to generate data for
visualizing the convergence of network states. The process involves two main
phases for each environmental condition and threshold:

1.  **Reference Attractor Generation:** A single simulation is run to find a
    stable fixed-point attractor, which serves as a reference state.
2.  **Population Convergence Simulation:** A large population of initial states
    is simulated asynchronously. At each epoch, the Hamming distance of every
    state in the population to the reference attractor is calculated.

The output is a set of compressed NumPy files (`.npz`) containing the full
trajectories of Hamming distances, suitable for generating violin plots.
'''

import os
import time
import numpy as np
import pandas as pd
from numba import njit
from pathlib import Path
from tqdm import tqdm

# =============================================================================
# SETUP: SIMULATION PARAMETERS AND FILE PATHS
# =============================================================================
RANDOM_SEED = 42
np.random.seed(RANDOM_SEED)

# Simulation Parameters
REALIZATIONS = 5_000_000
NUM_EPOCHS_BASE = 20
EXTRA_EPOCHS_PER_ATTEMPT = 5
MAX_CONVERGENCE_ATTEMPTS = 10
THRESHOLDS_TO_RUN = [2, 3]
INITIAL_ACTIVE_FRACTION = 0.5

# File Paths
BASE_DATA_PATH = Path("data")
OUTPUT_DIR = Path("violin_plot_data")
OUTPUT_DIR.mkdir(exist_ok=True)

# =============================================================================
# DATA LOADING AND PREPARATION
# =============================================================================

def load_network_data():
    """Loads all necessary network data and defines node groups."""
    print("--- 1. Loading and preparing network data ---")
    edge_list = pd.read_csv(BASE_DATA_PATH / "giantC_edge_list.csv")
    nodes_id = pd.read_csv(BASE_DATA_PATH / "giantC_nodes_id.csv")
    
    unique_ids = nodes_id['ID'].values
    id_to_index = {node_id: idx for idx, node_id in enumerate(unique_ids)}
    name_to_id = dict(zip(nodes_id['Gene'], nodes_id['ID']))
    N = len(unique_ids)
    
    # Filter edges to ensure all nodes exist in the mapping
    valid_edges = edge_list[edge_list['Node 1'].isin(id_to_index) & edge_list['Node 2'].isin(id_to_index)]
    src = np.array([id_to_index[i] for i in valid_edges['Node 1']], dtype=np.int32)
    dst = np.array([id_to_index[j] for i in valid_edges['Node 2']], dtype=np.int32)
    wts = valid_edges['Weight'].astype(np.int8).to_numpy()

    return N, id_to_index, name_to_id, src, dst, wts

def prepare_csr_representation(N, src, dst):
    """Prepares the network in a CSR-like format for Numba."""
    k_in = np.zeros(N, dtype=np.int32)
    for target_node in dst:
        k_in[target_node] += 1
        
    # An outgoing CSR is needed for efficient influence updates
    out_degree = np.zeros(N, dtype=np.int32)
    for source_node in src:
        out_degree[source_node] += 1
        
    row_ptr = np.concatenate((np.array([0]), np.cumsum(out_degree))).astype(np.int32)
    
    # Sort edges by source node to build the CSR structure correctly
    sort_indices = np.argsort(src)
    edges_dst = dst[sort_indices]
    edges_wt = wts[sort_indices]
    
    return k_in, row_ptr, edges_dst, edges_wt

# =============================================================================
# HIGH-PERFORMANCE SIMULATION FUNCTIONS (NUMBA-JIT)
# =============================================================================

@njit
def _run_single_trajectory(state, k_in, threshold, steps, src, dst, wts, out_row_ptr, out_edges_dst, out_edges_wt):
    """Core engine to simulate a single state vector for a number of steps."""
    N = state.shape[0]
    influence = np.zeros(N, dtype=np.int32)
    # Initial influence calculation
    for i in range(len(src)):
        if state[src[i]] == 1:
            influence[dst[i]] += wts[i]
            
    for _ in range(steps):
        j = np.random.randint(0, N)
        if k_in[j] > 0:
            new_state = 1 if influence[j] >= threshold else 0
            if new_state != state[j]:
                delta = new_state - state[j]
                state[j] = new_state
                # Update influence of neighbors efficiently using outgoing CSR
                start, end = out_row_ptr[j], out_row_ptr[j+1]
                for k in range(start, end):
                    neighbor = out_edges_dst[k]
                    weight = out_edges_wt[k]
                    influence[neighbor] += delta * weight
    return state

@njit
def _is_fixed_point(state, k_in, threshold, src, dst, wts):
    """Tests if a state is a fixed point under synchronous updates."""
    N = state.shape[0]
    influence = np.zeros(N, dtype=np.int32)
    for i in range(len(src)):
        if state[src[i]] == 1:
            influence[dst[i]] += wts[i]
    
    next_state = state.copy()
    for i in range(N):
        if k_in[i] > 0:
            next_state[i] = 1 if influence[i] >= threshold else 0
            
    return np.all(next_state == state)

@njit(parallel=True)
def calculate_hamming_distances(population, reference):
    """Calculates Hamming distances for a population against a reference state."""
    num_realizations, N = population.shape
    distances = np.zeros(num_realizations, dtype=np.int32)
    for r in range(num_realizations):
        distances[r] = np.sum(population[r] != reference)
    return distances

# =============================================================================
# MAIN SIMULATION WORKFLOW
# =============================================================================

def find_reference_attractor(env_params, N, k_in, threshold, non_fixed_indices, *csr_args):
    """Finds a stable reference attractor for a given environment."""
    print("    Phase 1: Generating reference attractor...")
    initial_state = np.zeros(N, dtype=np.int8)
    nodes_to_activate = np.random.choice(non_fixed_indices, size=int(len(non_fixed_indices) * INITIAL_ACTIVE_FRACTION), replace=False)
    initial_state[nodes_to_activate] = 1
    initial_state[env_params['active_indices']] = 1

    current_state = initial_state
    for attempt in range(MAX_CONVERGENCE_ATTEMPTS):
        steps = (NUM_EPOCHS_BASE * N) if attempt == 0 else (EXTRA_EPOCHS_PER_ATTEMPT * N)
        current_state = _run_single_trajectory(current_state, k_in, threshold, steps, *csr_args)
        if _is_fixed_point(current_state, k_in, threshold, src, dst, wts):
            print("      -> Reference attractor found to be a stable fixed point.")
            return current_state
            
    print("      -> [WARNING] Reference attractor did not stabilize into a fixed point.")
    return current_state

def simulate_population_convergence(ref_attractor, env_params, N, k_in, threshold, non_fixed_indices, *csr_args):
    """Simulates the population and records Hamming distance trajectories."""
    print("    Phase 2: Simulating population convergence...")
    
    # Create initial population
    population = np.zeros((REALIZATIONS, N), dtype=np.int8)
    for r in tqdm(range(REALIZATIONS), desc="      Initializing population"):
        nodes_to_activate = np.random.choice(non_fixed_indices, size=int(len(non_fixed_indices) * INITIAL_ACTIVE_FRACTION), replace=False)
        population[r, nodes_to_activate] = 1
        population[r, env_params['active_indices']] = 1
    
    # Simulate epoch by epoch
    hamming_trajectories = [calculate_hamming_distances(population, ref_attractor)]
    pbar = tqdm(total=NUM_EPOCHS_BASE, desc="      Simulating epochs")
    for epoch in range(NUM_EPOCHS_BASE):
        for r in range(REALIZATIONS):
            population[r] = _run_single_trajectory(population[r], k_in, threshold, N, *csr_args)
        hamming_trajectories.append(calculate_hamming_distances(population, ref_attractor))
        pbar.update(1)
    pbar.close()
    
    return np.array(hamming_trajectories)

# =============================================================================
# SCRIPT EXECUTION
# =============================================================================

if __name__ == "__main__":
    total_start_time = time.time()
    
    # 1. Load and prepare all data structures
    N_nodes, id_to_idx, name_to_id_map, src, dst, wts = load_network_data()
    k_in, out_row_ptr, out_edges_dst, out_edges_wt = prepare_csr_representation(N_nodes, src, dst)
    csr_arguments = (src, dst, wts, out_row_ptr, out_edges_dst, out_edges_wt)
    
    # Define environments
    def _get_indices(gene_list):
        return [id_to_idx[name_to_id_map[name]] for name in gene_list if name in name_to_id_map]
    
    environments = {
        "YPD": {'active_indices': np.array(_get_indices(["YDL056W", "YER111C"]))}, # Add full lists
        "Starvation": {'active_indices': np.array(_get_indices(["YMR037C", "YGL073W"]))} # Add full lists
    }
    non_fixed = np.where(k_in > 0)[0] # Simplified logic for non-fixed nodes

    # 2. Run simulation for each condition
    for t in THRESHOLDS_TO_RUN:
        print(f"\n--- PROCESSING FOR THRESHOLD t={t} ---")
        trajectory_data = {}
        for env_name, params in environments.items():
            print(f"  Environment: {env_name}")
            ref_state = find_reference_attractor(params, N_nodes, k_in, t, non_fixed, *csr_arguments)
            trajectories = simulate_population_convergence(ref_state, params, N_nodes, k_in, t, non_fixed, *csr_arguments)
            trajectory_data[f"{env_name.lower()}_trajectories"] = trajectories
        
        # 3. Save results for the current threshold
        output_path = OUTPUT_DIR / f"convergence_trajectories_t{t}.npz"
        np.savez_compressed(output_path, **trajectory_data)
        print(f"  -> Results for t={t} saved to '{output_path}'")

    print(f"\nTotal execution time: {(time.time() - total_start_time) / 60:.2f} minutes.")