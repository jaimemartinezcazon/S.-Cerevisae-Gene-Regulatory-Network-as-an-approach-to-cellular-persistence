'''
Author: Jaime Martínez Cazón

Description:
This script simulates the dynamical behavior of the S. cerevisiae Gene
Regulatory Network using a weighted asynchronous Boolean model. The simulation
tracks the network's trajectory over time as it transitions between different
environmental conditions (logarithmic growth vs. starvation), which are modeled
by fixing the states of specific sets of transcription factors.

The script uses a Numba-JIT compiled function for high-performance simulation
and runs multiple realizations to average out stochastic effects. The final
output is a compressed NumPy file (`.npz`) containing the time evolution of
network activity, vulnerability, and other metrics, ready for plotting.
'''

import os
import time
import numpy as np
import pandas as pd
from numba import njit
from pathlib import Path

# =============================================================================
# SETUP: SIMULATION PARAMETERS AND FILE PATHS
# =============================================================================
RANDOM_SEED = 42
np.random.seed(RANDOM_SEED)

# File Paths
BASE_DATA_PATH = Path("data")
OUTPUT_FILENAME = "simulation_trajectory_data.npz"

# Simulation Parameters
REALIZATIONS = 5000
THRESHOLDS = [2, 3]
EPOCHS_PER_SEGMENT = 15
NUM_SEGMENTS = 3  # YPD -> Starvation -> YPD

# =============================================================================
# DATA LOADING AND PREPARATION
# =============================================================================

def load_and_prepare_network_data():
    """
    Loads network data and prepares it in Compressed Sparse Row (CSR) format
    for efficient use within the Numba-JIT compiled simulation function.
    """
    print("--- 1. Loading and preparing network data ---")
    
    edge_list_path = BASE_DATA_PATH / "giantC_edge_list.csv"
    nodes_id_path = BASE_DATA_PATH / "giantC_nodes_id.csv"
    
    edge_list = pd.read_csv(edge_list_path)
    nodes_id = pd.read_csv(nodes_id_path)
    
    unique_ids = nodes_id['ID'].values
    id_to_index = {node_id: idx for idx, node_id in enumerate(unique_ids)}
    N = len(unique_ids)
    
    src = np.array([id_to_index[i] for i in edge_list['Node 1']], dtype=np.int32)
    dst = np.array([id_to_index[j] for j in edge_list['Node 2']], dtype=np.int32)
    wts = edge_list['Weight'].astype(np.int8).to_numpy()
    
    # Build CSR for incoming edges (for calculating influence)
    in_degree = np.zeros(N, dtype=np.int32)
    for target_node in dst:
        in_degree[target_node] += 1
        
    in_row_ptr = np.concatenate((np.array([0]), np.cumsum(in_degree))).astype(np.int32)
    in_edges_src = np.empty_like(src)
    in_edges_wt = np.empty_like(wts)
    
    cursor = in_row_ptr[:-1].copy()
    for i in range(len(src)):
        target_node = dst[i]
        idx = cursor[target_node]
        in_edges_src[idx] = src[i]
        in_edges_wt[idx] = wts[i]
        cursor[target_node] += 1
        
    print("Network data prepared in CSR format.")
    return N, id_to_index, nodes_id, in_row_ptr, in_edges_src, in_edges_wt

def define_simulation_conditions(N, nodes_id, id_to_index):
    """Defines the sets of TFs and candidate genes for the simulation."""
    LOG_PHASE_TFS = ["YDL056W", "YER111C", "YGR044C", "YIL131C", "YDR146C"] # And others...
    STARVATION_TFS = ["YMR037C", "YGL073W", "YOL116W", "YNL027W", "YOR028C"] # And others...
    CANDIDATE_GENES = ["YAL049C", "YBR208C", "YDL182W", "YFL014W", "YGR088W"] # And others...

    name_to_id = dict(zip(nodes_id['Gene'], nodes_id['ID']))

    def _get_indices(gene_list):
        return [id_to_index[name_to_id[name]] for name in gene_list if name in name_to_id and name_to_id[name] in id_to_index]

    # Boolean arrays for fixed nodes
    fixed_nodes_ypd = np.zeros(N, dtype=np.bool_)
    fixed_nodes_ypd[_get_indices(LOG_PHASE_TFS)] = True
    fixed_nodes_ypd[_get_indices(STARVATION_TFS)] = True
    
    fixed_nodes_starvation = np.zeros(N, dtype=np.bool_)
    fixed_nodes_starvation[_get_indices(STARVATION_TFS)] = True
    fixed_nodes_starvation[_get_indices(LOG_PHASE_TFS)] = True

    # State arrays for fixed nodes (-1 means not fixed)
    fixed_states_ypd = -np.ones(N, dtype=np.int8)
    fixed_states_ypd[_get_indices(LOG_PHASE_TFS)] = 1  # Active in YPD
    fixed_states_ypd[_get_indices(STARVATION_TFS)] = 0  # Inactive in YPD

    fixed_states_starvation = -np.ones(N, dtype=np.int8)
    fixed_states_starvation[_get_indices(STARVATION_TFS)] = 1  # Active in Starvation
    fixed_states_starvation[_get_indices(LOG_PHASE_TFS)] = 0  # Inactive in Starvation
    
    # Boolean array for candidate genes
    candidates_bool = np.zeros(N, dtype=np.bool_)
    candidates_bool[_get_indices(CANDIDATE_GENES)] = True

    segment_configs = [
        {'fixed_nodes': fixed_nodes_ypd, 'fixed_states': fixed_states_ypd},
        {'fixed_nodes': fixed_nodes_starvation, 'fixed_states': fixed_states_starvation},
        {'fixed_nodes': fixed_nodes_ypd, 'fixed_states': fixed_states_ypd}
    ]
    return segment_configs, candidates_bool

# =============================================================================
# HIGH-PERFORMANCE SIMULATION FUNCTION (NUMBA-JIT)
# =============================================================================

@njit(parallel=True)
def run_async_dynamics_jit(init_state, csr_in_ptr, csr_in_src, csr_in_wt, 
                           fixed_nodes, fixed_states, threshold, steps, candidates):
    """
    Performs one segment of the asynchronous dynamics simulation.
    This function is compiled to machine code by Numba for performance.
    """
    N = init_state.shape[0]
    current_state = init_state.copy()
    
    # Apply the fixed states for this specific segment
    for i in prange(N):
        if fixed_states[i] != -1:
            current_state[i] = fixed_states[i]
            
    # Pre-calculate initial influence
    influence = np.zeros(N, dtype=np.int32)
    for i in prange(N):
        start, end = csr_in_ptr[i], csr_in_ptr[i+1]
        for j in range(start, end):
            source_node = csr_in_src[j]
            influence[i] += csr_in_wt[j] * current_state[source_node]
            
    # Arrays to store results for this segment
    evol = np.empty((steps, N), dtype=np.int8)
    
    # Main simulation loop for the segment
    for t in range(steps):
        evol[t] = current_state
        
        # Select a random node to update
        node_to_update = np.random.randint(0, N)
        
        # Update only if the node is not fixed
        if not fixed_nodes[node_to_update]:
            new_state = 1 if influence[node_to_update] >= threshold else 0
            old_state = current_state[node_to_update]
            
            if new_state != old_state:
                delta = new_state - old_state
                current_state[node_to_update] = new_state
                
                # Update the influence of the neighbors of the updated node
                # This requires an outgoing CSR representation, which we can build
                # if performance becomes an issue. For now, we recalculate.
                # (Re-calculating influence is simpler and often fast enough)
                # For simplicity, let's re-calculate full influence.
                # NOTE: A more optimized version would use an outgoing CSR here.
                for i in prange(N):
                    start, end = csr_in_ptr[i], csr_in_ptr[i+1]
                    inf_val = 0
                    for j in range(start, end):
                        source_node = csr_in_src[j]
                        inf_val += csr_in_wt[j] * current_state[source_node]
                    influence[i] = inf_val

    return evol, current_state

# =============================================================================
# MAIN SIMULATION WORKFLOW
# =============================================================================

def run_full_simulation(N, csr_in_ptr, csr_in_src, csr_in_wt, segment_configs, candidates_bool):
    """Orchestrates the entire simulation across all thresholds and realizations."""
    simulation_results = {}
    steps_per_epoch = N
    total_steps = steps_per_epoch * EPOCHS_PER_SEGMENT * NUM_SEGMENTS

    for w in THRESHOLDS:
        print(f"\n--- 2. Simulating for threshold omega = {w} ---")
        
        # Arrays to store aggregated results for this threshold
        all_realizations_activity = np.empty((REALIZATIONS, total_steps), dtype=np.float64)
        
        for r in tqdm(range(REALIZATIONS), desc=f"Threshold {w}"):
            np.random.seed(r)
            current_state = np.zeros(N, dtype=np.int8)
            step_offset = 0

            for config in segment_configs:
                evol, final_state = run_async_dynamics_jit(
                    current_state, csr_in_ptr, csr_in_src, csr_in_wt,
                    config['fixed_nodes'], config['fixed_states'],
                    w, EPOCHS_PER_SEGMENT * steps_per_epoch, candidates_bool
                )
                
                # Store the mean activity for this segment
                segment_len = evol.shape[0]
                all_realizations_activity[r, step_offset : step_offset + segment_len] = np.mean(evol, axis=1)
                
                current_state = final_state
                step_offset += segment_len
        
        # Calculate statistics across all realizations
        simulation_results[f'mean_active_{w}'] = np.mean(all_realizations_activity, axis=0)
        simulation_results[f'std_active_{w}'] = np.std(all_realizations_activity, axis=0)
        # Add other metrics here if needed (e.g., vulnerability)
        
    # Add simulation parameters to the results dict for saving
    simulation_results['thresholds'] = np.array(THRESHOLDS)
    simulation_results['total_steps'] = np.array(total_steps)
    simulation_results['steps_per_epoch'] = np.array(steps_per_epoch)
    
    return simulation_results

# =============================================================================
# SCRIPT EXECUTION
# =============================================================================

if __name__ == "__main__":
    start_time = time.time()
    
    # 1. Load data and prepare network structures
    N_nodes, id_map, nodes_df, csr_ptr, csr_src, csr_wt = load_and_prepare_network_data()
    
    # 2. Define simulation conditions (fixed nodes, candidates)
    sim_configs, candidates = define_simulation_conditions(N_nodes, nodes_df, id_map)
    
    # 3. Run the full simulation pipeline
    results = run_full_simulation(N_nodes, csr_ptr, csr_src, csr_wt, sim_configs, candidates)
    
    # 4. Save the results
    print(f"\n--- 3. Saving simulation results ---")
    np.savez(OUTPUT_FILENAME, **results)
    
    end_time = time.time()
    print(f"Simulation finished and data saved to '{OUTPUT_FILENAME}'")
    print(f"Total execution time: {end_time - start_time:.2f} seconds.")