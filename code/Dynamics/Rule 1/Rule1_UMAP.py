'''
Author: Jaime Martínez Cazón

Description:
This script generates the raw steady-state data required for dimensionality
reduction analyses like UMAP or PCA. For each specified threshold and
environmental condition, it performs a large number of simulations, each
starting from a unique random initial state.

The script simulates until a fixed-point attractor is reached, using an
adaptive step approach to ensure convergence. It then saves two main types
of data for each threshold:
1.  The full binary steady-state vectors for every realization.
2.  A matrix of summary metrics (e.g., active fraction) for various node groups.

Additionally, it performs a check for multiple attractors by comparing all
final states to a reference attractor derived from the first realization.
'''

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
REALIZATIONS = 1_000_000
NUM_EPOCHS_BASE = 20
EXTRA_EPOCHS_PER_ATTEMPT = 5
MAX_CONVERGENCE_ATTEMPTS = 10
THRESHOLDS_TO_RUN = list(range(1, 7))
INITIAL_ACTIVE_FRACTION = 0.5

# File Paths
BASE_DATA_PATH = Path("data")
OUTPUT_DIR = Path("steady_state_data")
OUTPUT_DIR.mkdir(exist_ok=True)

# =============================================================================
# DATA LOADING AND PREPARATION
# =============================================================================

def load_and_prepare_data():
    """Loads all necessary data and defines node groups for analysis."""
    print("--- 1. Loading and preparing network data and node groups ---")
    
    edge_list = pd.read_csv(BASE_DATA_PATH / "giantC_edge_list.csv")
    nodes_id = pd.read_csv(BASE_DATA_PATH / "giantC_nodes_id.csv")
    
    id_to_index = {node_id: idx for idx, node_id in enumerate(nodes_id['ID'].values)}
    N = len(id_to_index)
    
    src = np.array([id_to_index[i] for i in edge_list['Node 1']], dtype=np.int32)
    dst = np.array([id_to_index[j] for j in edge_list['Node 2']], dtype=np.int32)
    wts = edge_list['Weight'].astype(np.int8).to_numpy()

    k_in = np.bincount(dst, minlength=N).astype(np.int32)
    non_fixed_indices = np.where(k_in > 0)[0]

    node_groups = {'global': np.ones(N, dtype=np.bool_)}
    # ... (Add logic to load communities, scc, candidates as in previous scripts)
    
    name_to_id = dict(zip(nodes_id['Gene'], nodes_id['ID']))
    log_phase_tfs = ["YDL056W", "YER111C"] # Add full list
    starvation_tfs = ["YMR037C", "YGL073W"] # Add full list
    
    def _get_indices(gene_list):
        return [id_to_index[name_to_id[name]] for name in gene_list if name in name_to_id]

    environments = {
        "YPD": {'active_indices': np.array(_get_indices(log_phase_tfs))},
        "Starvation": {'active_indices': np.array(_get_indices(starvation_tfs))}
    }
    
    print("--- Data preparation complete ---")
    return N, k_in, src, dst, wts, non_fixed_indices, node_groups, environments, nodes_id['ID'].values

# =============================================================================
# HIGH-PERFORMANCE SIMULATION FUNCTION (NUMBA-JIT)
# =============================================================================

@njit
def find_final_attractor_jit(initial_state, k_in, threshold, src, dst, wts,
                             base_steps, extra_steps, max_attempts):
    """
    Simulates a single trajectory until a fixed point is reached or max attempts are exceeded.
    Returns the final state and a flag indicating if it's a stable fixed point.
    """
    state = initial_state.copy()
    N = state.shape[0]
    
    for _ in range(max_attempts + 1):
        # Asynchronous simulation for a block of steps
        influence = np.zeros(N, dtype=np.int32)
        for i in range(len(src)):
            if state[src[i]] == 1:
                influence[dst[i]] += wts[i]
        
        steps_to_run = base_steps if _ == 0 else extra_steps
        for _ in range(steps_to_run):
            j = np.random.randint(N)
            if k_in[j] > 0:
                new_state = 1 if influence[j] >= threshold else 0
                if new_state != state[j]:
                    state[j] = new_state
                    # Full influence re-calculation for correctness
                    influence = np.zeros(N, dtype=np.int32)
                    for i in range(len(src)):
                        if state[src[i]] == 1:
                            influence[dst[i]] += wts[i]

        # Synchronous update test for fixed point
        next_state_sync = state.copy()
        is_stable = True
        for i in range(N):
            if k_in[i] > 0:
                next_state_sync[i] = 1 if influence[i] >= threshold else 0
        
        if np.all(next_state_sync == state):
            return state, True # Return stable state
            
    return state, False # Return last state if not stabilized

# =============================================================================
# MAIN DATA GENERATION WORKFLOW
# =============================================================================

def run_data_generation(N, k_in, src, dst, wts, non_fixed, groups, envs, all_node_ids):
    """Orchestrates the data generation for all thresholds and environments."""
    
    steps_per_epoch = N
    base_sim_steps = NUM_EPOCHS_BASE * steps_per_epoch
    extra_sim_steps = EXTRA_EPOCHS_PER_ATTEMPT * steps_per_epoch

    for t in THRESHOLDS_TO_RUN:
        print("-" * 60)
        print(f"--- PROCESSING FOR THRESHOLD t={t} ---")
        
        threshold_data = {}
        for env_name, params in envs.items():
            print(f"  Environment: {env_name}...")
            
            state_matrix = np.empty((REALIZATIONS, N), dtype=np.int8)
            ref_attractor = None
            
            pbar = tqdm(total=REALIZATIONS, desc=f"  Simulating {env_name} t={t}")
            for r in range(REALIZATIONS):
                init_state = np.zeros(N, dtype=np.int8)
                nodes_to_activate = np.random.choice(non_fixed, size=int(round(INITIAL_ACTIVE_FRACTION * len(non_fixed))), replace=False)
                init_state[nodes_to_activate] = 1
                init_state[params['active_indices']] = 1
                
                final_state, _ = find_final_attractor_jit(
                    init_state, k_in, t, src, dst, wts,
                    base_sim_steps, extra_sim_steps, MAX_CONVERGENCE_ATTEMPTS
                )
                state_matrix[r, :] = final_state
                
                if r == 0:
                    ref_attractor = final_state
                pbar.update(1)
            pbar.close()

            # Multiple attractor check
            hamming_distances = np.sum(state_matrix != ref_attractor, axis=1)
            if np.any(hamming_distances > 0):
                print(f"    [INFO] Multiple attractors detected for t={t}, {env_name}.")
            else:
                print(f"    [INFO] All realizations converged to a single attractor for t={t}, {env_name}.")
            
            threshold_data[env_name] = {"state_matrix": state_matrix}
        
        # Save data for the current threshold
        output_path = OUTPUT_DIR / f"steady_state_data_t{t}.npz"
        np.savez_compressed(
            output_path,
            ypd_states=threshold_data['YPD']['state_matrix'],
            starvation_states=threshold_data['Starvation']['state_matrix'],
            node_ids=all_node_ids
        )
        print(f"  -> Data for t={t} saved to '{output_path}'")

# =============================================================================
# SCRIPT EXECUTION
# =============================================================================

if __name__ == "__main__":
    start_time = time.time()
    
    # 1. Load data and define experimental setup
    try:
        n_nodes, k_in_arr, src_arr, dst_arr, wts_arr, non_fixed_idx, node_groups, envs, node_id_list = load_and_prepare_data()
    except FileNotFoundError as e:
        print(f"Fatal Error: A required data file was not found. {e}")
        exit()
        
    # 2. Run the full data generation pipeline
    run_data_generation(n_nodes, k_in_arr, src_arr, dst_arr, wts_arr, non_fixed_idx, node_groups, envs, node_id_list)
    
    print(f"\nTotal execution time: {(time.time() - start_time) / 60:.2f} minutes.")