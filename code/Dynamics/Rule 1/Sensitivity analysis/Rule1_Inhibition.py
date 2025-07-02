'''
Author: Jaime Martínez Cazón

Description:
This script performs a sensitivity analysis on the Gene Regulatory Network's
steady-state behavior by varying the percentage of inhibitory connections.
For each realization, the script ignores the original network weights and
randomly reassigns a specified fraction of edges to be inhibitory (-1), with
the rest being activating (+1).

It systematically explores a parameter space defined by the activation
threshold and the inhibition percentage, running multiple simulations to find
the average final state for various node groups. This analysis tests the
robustness of the network's dynamics to changes in its regulatory logic.
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
NUM_EPOCHS_BASE = 20
EXTRA_EPOCHS_PER_ATTEMPT = 5
MAX_CONVERGENCE_ATTEMPTS = 10
THRESHOLDS = [1, 2, 3, 4, 5, 6]
INHIBITION_PERCS = np.arange(0, 101, 5) / 100.0
INITIAL_ACTIVE_FRACTION = 0.5

# File Paths
BASE_DATA_PATH = Path("data")
OUTPUT_DIR = Path("inhibition_sensitivity_output")
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
    NUM_EDGES = len(edge_list)
    
    src = np.array([id_to_index[i] for i in edge_list['Node 1']], dtype=np.int32)
    dst = np.array([id_to_index[j] for j in edge_list['Node 2']], dtype=np.int32)
    
    k_in = np.bincount(dst, minlength=N).astype(np.int32)
    non_input_indices = np.where(k_in > 0)[0]

    node_groups = {'global': np.ones(N, dtype=np.bool_)}
    # ... (Add loading for candidates, communities, SCC as before)

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
    return N, NUM_EDGES, k_in, src, dst, non_input_indices, node_groups, environments

# =============================================================================
# HIGH-PERFORMANCE SIMULATION FUNCTION (NUMBA-JIT)
# =============================================================================

@njit
def run_single_realization_jit(initial_state, k_in, threshold, p_inhib,
                               num_edges, src, dst, base_steps, extra_steps, max_attempts, group_masks):
    """
    Core JIT-compiled engine for one full realization.
    It randomizes weights, simulates until convergence, and returns final fractions.
    """
    # 1. Randomly assign weights for this realization
    wts = np.ones(num_edges, dtype=np.int8)
    num_inhibitory = int(round(p_inhib * num_edges))
    if num_inhibitory > 0:
        inhib_indices = np.random.choice(np.arange(num_edges), size=num_inhibitory, replace=False)
        wts[inhib_indices] = -1

    # 2. Adaptive simulation loop
    state = initial_state.copy()
    N = state.shape[0]
    is_stable = False
    
    for attempt in range(max_attempts + 1):
        steps_to_run = base_steps if attempt == 0 else extra_steps
        
        influence = np.zeros(N, dtype=np.int32)
        for i in range(num_edges):
            if state[src[i]] == 1:
                influence[dst[i]] += wts[i]
                
        for _ in range(steps_to_run):
            j = np.random.randint(N)
            if k_in[j] > 0:
                new_state = 1 if influence[j] >= threshold else 0
                if new_state != state[j]:
                    # Recalculate full influence after a state change for simplicity and correctness
                    state[j] = new_state
                    influence = np.zeros(N, dtype=np.int32)
                    for i in range(num_edges):
                        if state[src[i]] == 1:
                            influence[dst[i]] += wts[i]

        # Test for stability (fixed point)
        next_state_sync = state.copy()
        for i in range(N):
            if k_in[i] > 0:
                next_state_sync[i] = 1 if influence[i] >= threshold else 0
        
        if np.all(next_state_sync == state):
            is_stable = True
            break
    
    # 3. Calculate final fractions
    num_groups = len(group_masks)
    final_fractions = np.empty(num_groups * 3, dtype=np.float64)
    for i in range(num_groups):
        mask = group_masks[i]
        total_in_group = np.sum(mask)
        if total_in_group > 0:
            active = np.sum(state[mask]) / total_in_group
            # Vulnerability calculations can be added here if needed
            final_fractions[i*3] = active
            final_fractions[i*3 + 1] = 0.0 # Placeholder for On-Vuln
            final_fractions[i*3 + 2] = 0.0 # Placeholder for Off-Vuln
        else:
            final_fractions[i*3 : i*3 + 3] = 0.0

    return final_fractions, is_stable

# =============================================================================
# MAIN SIMULATION WORKFLOW
# =============================================================================

def run_inhibition_analysis(N, num_edges, k_in, src, dst, non_input, groups, envs):
    """Orchestrates the full simulation to generate the inhibition sensitivity data."""
    
    group_masks_tuple = tuple(groups.values())
    steps_per_epoch = N
    base_sim_steps = NUM_EPOCHS_BASE * steps_per_epoch
    extra_sim_steps = EXTRA_EPOCHS_PER_ATTEMPT * steps_per_epoch
    
    for env_name, params in envs.items():
        print(f"\n--- 2. Starting analysis for environment: {env_name} ---")
        
        results_mean = np.zeros((len(THRESHOLDS), len(INHIBITION_PERCS), len(groups) * 3))
        results_std = np.zeros_like(results_mean)
        
        pbar_outer = tqdm(total=len(THRESHOLDS) * len(INHIBITION_PERCS), desc=f"Env: {env_name}")

        for i_t, t in enumerate(THRESHOLDS):
            for i_p, p_inhib in enumerate(INHIBITION_PERCS):
                realization_results = np.empty((REALIZATIONS, len(groups) * 3))
                unstable_count = 0
                
                for r in range(REALIZATIONS):
                    # Create initial state (same for all realizations in this inner loop)
                    init_state = np.zeros(N, dtype=np.int8)
                    nodes_to_activate = np.random.choice(non_input, size=int(round(INITIAL_ACTIVE_FRACTION * len(non_input))), replace=False)
                    init_state[nodes_to_activate] = 1
                    init_state[params['active_indices']] = 1
                    
                    fractions, is_stable = run_single_realization_jit(
                        init_state, k_in, t, p_inhib, num_edges, src, dst,
                        base_sim_steps, extra_sim_steps, MAX_CONVERGENCE_ATTEMPTS, group_masks_tuple
                    )
                    realization_results[r, :] = fractions
                    if not is_stable: unstable_count += 1
                
                results_mean[i_t, i_p, :] = np.mean(realization_results, axis=0)
                results_std[i_t, i_p, :] = np.std(realization_results, axis=0)
                pbar_outer.update(1)
        
        pbar_outer.close()
        
        # Save results for this environment
        output_path = OUTPUT_DIR / f"{env_name}_inhibition_sensitivity_data.npz"
        np.savez_compressed(
            output_path, thresholds=np.array(THRESHOLDS), inhibition_percs=INHIBITION_PERCS,
            results_mean=results_mean, results_std=results_std, group_names=np.array(list(groups.keys()))
        )
        print(f"\nResults for {env_name} saved to '{output_path}'")

# =============================================================================
# SCRIPT EXECUTION
# =============================================================================

if __name__ == "__main__":
    start_time = time.time()
    
    # 1. Load data and define setup
    N_nodes, n_edges, k_in_arr, src_arr, dst_arr, non_input_idx, node_groups_map, envs_map = load_and_prepare_data()
    
    # 2. Run the inhibition sensitivity analysis
    run_inhibition_analysis(N_nodes, n_edges, k_in_arr, src_arr, dst_arr, non_input_idx, node_groups_map, envs_map)
    
    print(f"\nTotal execution time: {(time.time() - start_time) / 60:.2f} minutes.")
