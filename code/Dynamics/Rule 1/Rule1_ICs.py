'''
Author: Jaime Martínez Cazón

Description:
This script generates a phase map of the Gene Regulatory Network's steady-state
behavior. It systematically explores a parameter space defined by two key
variables: the initial fraction of active nodes and the activation threshold.

For each combination of parameters and for different environmental conditions
(YPD vs. Starvation), the script runs multiple simulation realizations to find
the average final state. The output is a set of compressed NumPy files (`.npz`)
containing the mean and standard deviation of final state fractions for
various node groups, suitable for creating bifurcation diagrams or heatmaps.
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
REALIZATIONS = 5000
NUM_EPOCHS = 20
THRESHOLDS = [1, 2, 3, 4, 5]
INITIAL_COND_PERCS = np.arange(0, 101, 5) / 100.0

# File Paths
BASE_DATA_PATH = Path("data")
OUTPUT_DIR = Path("phase_map_output")
OUTPUT_DIR.mkdir(exist_ok=True)

# =============================================================================
# DATA LOADING AND PREPARATION
# =============================================================================

def load_and_prepare_data():
    """Loads all necessary data and defines node groups for analysis."""
    print("--- 1. Loading and preparing network data and node groups ---")
    
    # Load core network files
    edge_list = pd.read_csv(BASE_DATA_PATH / "giantC_edge_list.csv")
    nodes_id = pd.read_csv(BASE_DATA_PATH / "giantC_nodes_id.csv")
    
    id_to_index = {node_id: idx for idx, node_id in enumerate(nodes_id['ID'].values)}
    N = len(id_to_index)
    
    # Prepare network arrays
    src = np.array([id_to_index[i] for i in edge_list['Node 1']], dtype=np.int32)
    dst = np.array([id_to_index[j] for j in edge_list['Node 2']], dtype=np.int32)
    wts = edge_list['Weight'].astype(np.int8).to_numpy()
    
    k_in = np.bincount(dst, minlength=N).astype(np.int32)
    non_input_indices = np.where(k_in > 0)[0]

    # Define node groups
    node_groups = {'global': np.ones(N, dtype=np.bool_)}
    # ... (Add loading for candidates, communities, SCC as before)
    
    # Define environments
    name_to_id = dict(zip(nodes_id['Gene'], nodes_id['ID']))
    log_phase_tfs = ["YDL056W", "YER111C"] # Add full list
    starvation_tfs = ["YMR037C", "YGL073W"] # Add full list
    
    def _get_indices(gene_list):
        return [id_to_index[name_to_id[name]] for name in gene_list if name in name_to_id]

    environments = {
        "YPD": {'active_indices': np.array(_get_indices(log_phase_tfs), dtype=np.int32)},
        "Starvation": {'active_indices': np.array(_get_indices(starvation_tfs), dtype=np.int32)}
    }
    
    print("--- Data preparation complete ---")
    return N, k_in, src, dst, wts, non_input_indices, node_groups, environments

# =============================================================================
# HIGH-PERFORMANCE SIMULATION FUNCTION (NUMBA-JIT)
# =============================================================================

@njit
def run_simulation_engine_jit(initial_state, k_in, threshold, steps, src, dst, wts, group_masks):
    """
    Core JIT-compiled simulation engine. Simulates one realization and calculates
    final state fractions for all specified node groups.
    """
    state = initial_state.copy()
    N = state.shape[0]
    
    # Initial influence calculation
    influence = np.zeros(N, dtype=np.int32)
    for i in range(len(src)):
        if state[src[i]] == 1:
            influence[dst[i]] += wts[i]
            
    # Asynchronous dynamics loop
    for _ in range(steps):
        j = np.random.randint(N)
        if k_in[j] > 0:
            new_state = 1 if influence[j] >= threshold else 0
            if new_state != state[j]:
                # This simplified influence update is faster if k_out is small
                state[j] = new_state
                # A full re-calculation is needed for correctness after each update
                influence = np.zeros(N, dtype=np.int32)
                for i in range(len(src)):
                    if state[src[i]] == 1:
                        influence[dst[i]] += wts[i]
                        
    # Calculate final fractions for all groups
    num_groups = len(group_masks)
    final_fractions = np.empty(num_groups * 3, dtype=np.float64) # 3 metrics per group
    for i in range(num_groups):
        mask = group_masks[i]
        total_in_group = np.sum(mask)
        if total_in_group > 0:
            active_count, on_vuln, off_vuln = 0, 0, 0
            for node_idx in range(N):
                if mask[node_idx]:
                    active_count += state[node_idx]
                    if state[node_idx] == 0 and influence[node_idx] == threshold - 1:
                        on_vuln += 1
                    elif state[node_idx] == 1 and influence[node_idx] == threshold:
                        off_vuln += 1
            final_fractions[i*3] = active_count / total_in_group
            final_fractions[i*3 + 1] = on_vuln / total_in_group
            final_fractions[i*3 + 2] = off_vuln / total_in_group
        else:
            final_fractions[i*3 : i*3 + 3] = 0.0

    return final_fractions

# =============================================================================
# MAIN SIMULATION WORKFLOW
# =============================================================================

def run_phase_map_analysis(N, k_in, src, dst, wts, non_input, groups, envs):
    """Orchestrates the full simulation to generate the phase map data."""
    
    group_masks_tuple = tuple(groups.values())
    total_steps = NUM_EPOCHS * N
    
    for env_name, params in envs.items():
        print(f"\n--- 2. Starting analysis for environment: {env_name} ---")
        
        # Results arrays for this environment
        results_mean = np.zeros((len(THRESHOLDS), len(INITIAL_COND_PERCS), len(groups) * 3))
        results_std = np.zeros_like(results_mean)

        pbar = tqdm(total=len(THRESHOLDS) * len(INITIAL_COND_PERCS), desc=f"Env: {env_name}")
        
        for i_t, t in enumerate(THRESHOLDS):
            for i_p, p in enumerate(INITIAL_COND_PERCS):
                realization_results = np.empty((REALIZATIONS, len(groups) * 3))
                
                for r in range(REALIZATIONS):
                    # Create a new random initial state for each realization
                    init_state = np.zeros(N, dtype=np.int8)
                    num_to_activate = int(round(p * len(non_input)))
                    nodes_to_activate = np.random.choice(non_input, size=num_to_activate, replace=False)
                    init_state[nodes_to_activate] = 1
                    
                    # Enforce environmental conditions
                    init_state[params['active_indices']] = 1
                    
                    # Run simulation and store final fractions
                    realization_results[r, :] = run_simulation_engine_jit(
                        init_state, k_in, t, total_steps, src, dst, wts, group_masks_tuple
                    )
                
                results_mean[i_t, i_p, :] = np.mean(realization_results, axis=0)
                results_std[i_t, i_p, :] = np.std(realization_results, axis=0)
                pbar.update(1)
        
        pbar.close()
        
        # Save results for this environment
        output_path = OUTPUT_DIR / f"{env_name}_phase_map_data.npz"
        np.savez_compressed(
            output_path,
            thresholds=np.array(THRESHOLDS),
            initial_percs=INITIAL_COND_PERCS,
            results_mean=results_mean,
            results_std=results_std,
            group_names=np.array(list(groups.keys()))
        )
        print(f"\nResults for {env_name} saved to '{output_path}'")

# =============================================================================
# SCRIPT EXECUTION
# =============================================================================

if __name__ == "__main__":
    start_time = time.time()
    
    # 1. Load data and define experimental setup
    N_nodes, k_in_arr, src_arr, dst_arr, wts_arr, non_input_idx, node_groups_map, environments_map = load_and_prepare_data()
    
    # 2. Run the full phase map simulation
    run_phase_map_analysis(N_nodes, k_in_arr, src_arr, dst_arr, wts_arr, non_input_idx, node_groups_map, environments_map)
    
    print(f"\nTotal execution time: {(time.time() - start_time) / 60:.2f} minutes.")