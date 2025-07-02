'''
Author: Jaime Martínez Cazón

Description:
This script simulates the dynamical behavior of the Gene Regulatory Network using
a **proportional, node-specific threshold rule**. The activation threshold for
each node `j` is defined as `w * k_in(j)`, where `k_in(j)` is the node's
in-degree.

This experiment investigates how the network's trajectory changes as the
**proportionality constant `w`** is varied, while keeping the initial condition
fixed (50% active nodes). The results are saved to .npz files, one for each
tested value of `w`.
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
REALIZATIONS = 1000
NUM_EPOCHS = 15
INITIAL_ACTIVE_FRACTION = 0.5
PROPORTIONALITY_CONSTANTS_W = [0.3, 0.4, 0.5, 0.6, 0.7, 0.8]

# File Paths
BASE_DATA_PATH = Path("data")
OUTPUT_DIR = Path("proportional_threshold_vs_w_output")
OUTPUT_DIR.mkdir(exist_ok=True)

# =============================================================================
# DATA LOADING AND PREPARATION
# =============================================================================

def load_and_prepare_data():
    """Loads all necessary network data."""
    print("--- 1. Loading and preparing network data ---")
    
    edge_list = pd.read_csv(BASE_DATA_PATH / "giantC_edge_list.csv")
    nodes_id = pd.read_csv(BASE_DATA_PATH / "giantC_nodes_id.csv")
    
    id_to_index = {node_id: idx for idx, node_id in enumerate(nodes_id['ID'].values)}
    N = len(id_to_index)
    
    src = np.array([id_to_index[i] for i in edge_list['Node 1']], dtype=np.int32)
    dst = np.array([id_to_index[j] for j in edge_list['Node 2']], dtype=np.int32)
    wts = edge_list['Weight'].astype(np.int8).to_numpy()

    k_in = np.bincount(dst, minlength=N).astype(np.int32)
    
    # For this experiment, no nodes are fixed by environmental conditions
    fixed_nodes = np.zeros(N, dtype=np.bool_)
    
    print("--- Data preparation complete ---")
    return N, k_in, src, dst, wts, fixed_nodes

# =============================================================================
# HIGH-PERFORMANCE SIMULATION FUNCTION (NUMBA-JIT)
# =============================================================================

@njit
def run_proportional_threshold_jit(initial_state, k_in, threshold_array, steps, src, dst, wts, fixed_nodes):
    """
    Core JIT-compiled simulation engine using a pre-calculated array of
    proportional thresholds.
    """
    N = initial_state.shape[0]
    state = initial_state.copy()
    
    evol_active = np.empty(steps, dtype=np.float64)
    evol_on_vuln = np.empty(steps, dtype=np.float64)
    evol_off_vuln = np.empty(steps, dtype=np.float64)

    influence = np.zeros(N, dtype=np.int32)
    for i in range(len(src)):
        if state[src[i]] == 1:
            influence[dst[i]] += wts[i]
            
    for t in range(steps):
        active_count, on_vuln, off_vuln = 0, 0, 0
        for i in range(N):
            active_count += state[i]
            if state[i] == 0 and influence[i] == threshold_array[i] - 1: on_vuln += 1
            elif state[i] == 1 and influence[i] == threshold_array[i]: off_vuln += 1

        evol_active[t] = active_count / N
        evol_on_vuln[t] = on_vuln / N
        evol_off_vuln[t] = off_vuln / N

        j = np.random.randint(N)
        if not fixed_nodes[j] and k_in[j] > 0:
            new_state = 1 if influence[j] >= threshold_array[j] else 0
            if new_state != state[j]:
                state[j] = new_state
                influence = np.zeros(N, dtype=np.int32)
                for i in range(len(src)):
                    if state[src[i]] == 1:
                        influence[dst[i]] += wts[i]
                        
    return evol_active, evol_on_vuln, evol_off_vuln

# =============================================================================
# MAIN SIMULATION WORKFLOW
# =============================================================================

def run_analysis_for_w(w_factor, N, k_in, src, dst, wts, fixed_nodes):
    """
    Runs all realizations for a single proportionality constant 'w'.
    """
    print(f"  -> Running {REALIZATIONS} realizations for w = {w_factor:.2f}...")
    
    total_steps = NUM_EPOCHS * N
    
    # Pre-calculate the node-specific integer thresholds for this 'w'
    # Rule: threshold = max(1, floor(w * k_in))
    thresholds_per_node = np.maximum(1, np.floor(w_factor * k_in)).astype(np.int32)

    active_trajectories = np.empty((REALIZATIONS, total_steps))
    on_vuln_trajectories = np.empty((REALIZATIONS, total_steps))
    off_vuln_trajectories = np.empty((REALIZATIONS, total_steps))
    
    pbar = tqdm(range(REALIZATIONS), desc=f"w={w_factor:.2f}", leave=False)
    for r in pbar:
        init_state = np.random.choice([0, 1], size=N, p=[1 - INITIAL_ACTIVE_FRACTION, INITIAL_ACTIVE_FRACTION]).astype(np.int8)
        
        active_traj, on_vuln_traj, off_vuln_traj = run_proportional_threshold_jit(
            init_state, k_in, thresholds_per_node, total_steps, src, dst, wts, fixed_nodes
        )
        
        active_trajectories[r, :] = active_traj
        on_vuln_trajectories[r, :] = on_vuln_traj
        off_vuln_trajectories[r, :] = off_vuln_traj
        
    return {
        'm_act': np.mean(active_trajectories, axis=0), 's_act': np.std(active_trajectories, axis=0),
        'm_on': np.mean(on_vuln_trajectories, axis=0), 's_on': np.std(on_vuln_trajectories, axis=0),
        'm_off': np.mean(off_vuln_trajectories, axis=0), 's_off': np.std(off_vuln_trajectories, axis=0)
    }

# =============================================================================
# SCRIPT EXECUTION
# =============================================================================

if __name__ == "__main__":
    start_time = time.time()
    
    try:
        N_nodes, k_in_arr, src_arr, dst_arr, wts_arr, fixed_nodes_arr = load_and_prepare_data()
    except FileNotFoundError as e:
        print(f"Fatal Error: {e}")
        exit()

    print(f"\n--- 2. Starting simulation for different proportionality factors (w) ---")
    
    for w_val in PROPORTIONALITY_CONSTANTS_W:
        results = run_analysis_for_w(w_val, N_nodes, k_in_arr, src_arr, dst_arr, wts_arr, fixed_nodes_arr)
        
        output_filename = OUTPUT_DIR / f"proportional_threshold_pIC{INITIAL_ACTIVE_FRACTION:.1f}_w{w_val:.1f}.npz"
        np.savez(output_filename, w=w_val, **results)
        print(f"  -> Results saved to '{output_filename}'")
        
    print(f"\nTotal execution time: {(time.time() - start_time) / 60:.2f} minutes.")
