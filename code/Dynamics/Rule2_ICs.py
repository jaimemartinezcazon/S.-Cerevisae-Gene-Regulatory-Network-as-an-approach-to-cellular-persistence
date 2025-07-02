'''
Author: Jaime Martínez Cazón

Description:
This script simulates the dynamical behavior of the Gene Regulatory Network using
a modified asynchronous Boolean model. In this version, the activation
threshold for each node is not fixed but is **proportional** to its in-degree (k_in).
The update rule for a node `j` is: state becomes 1 if influence >= w * k_in(j),
where `w` is a fixed proportionality constant.

The script explores how the system's final state depends on the initial
fraction of active nodes for a given proportionality constant `w`.
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
PROPORTIONALITY_CONSTANT_W = 0.2  # The 'w' factor for the proportional threshold
REALIZATIONS = 1000
NUM_EPOCHS = 15
INITIAL_ACTIVE_FRACTIONS = [0.2, 0.4, 0.6, 0.8, 1.0] # p_one values

# File Paths
BASE_DATA_PATH = Path("data")
OUTPUT_DIR = Path("proportional_threshold_output")
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
    
    print("--- Data preparation complete ---")
    return N, k_in, src, dst, wts

# =============================================================================
# HIGH-PERFORMANCE SIMULATION FUNCTION (NUMBA-JIT)
# =============================================================================

@njit
def run_proportional_threshold_simulation_jit(initial_state, k_in, w, steps, src, dst, wts):
    """
    Core JIT-compiled simulation engine using a proportional threshold.
    Returns the full time evolution of activity and vulnerability fractions.
    """
    N = initial_state.shape[0]
    state = initial_state.copy()
    
    # Pre-calculate the integer threshold for each node
    # The threshold is ceil(w * k_in), but since influence is integer, it's equivalent
    # to influence >= w * k_in. For vulnerability check, we need integers.
    thresholds = (w * k_in).astype(np.int32) # For comparison. Actual threshold is float.
    
    # Store evolution of metrics
    evol_active = np.empty(steps, dtype=np.float64)
    evol_on_vuln = np.empty(steps, dtype=np.float64)
    evol_off_vuln = np.empty(steps, dtype=np.float64)

    # Initial influence calculation
    influence = np.zeros(N, dtype=np.int32)
    for i in range(len(src)):
        if state[src[i]] == 1:
            influence[dst[i]] += wts[i]
            
    # Asynchronous dynamics loop
    for t in range(steps):
        # Calculate and store metrics for the current state
        active_count, on_vuln, off_vuln = 0, 0, 0
        for i in range(N):
            active_count += state[i]
            # Vulnerability: influence is exactly one step away from crossing the threshold
            # Threshold for node i is w * k_in[i].
            # For integer influence, this means influence == floor(w * k_in[i]) for OFF
            # and influence == floor(w * k_in[i]) - 1 for ON
            # A simpler way is to use np.ceil
            node_threshold = np.ceil(w * k_in[i])
            if state[i] == 1 and influence[i] == node_threshold:
                off_vuln += 1
            elif state[i] == 0 and influence[i] == node_threshold - 1:
                on_vuln += 1

        evol_active[t] = active_count / N
        evol_on_vuln[t] = on_vuln / N
        evol_off_vuln[t] = off_vuln / N

        # Select a random node to update
        j = np.random.randint(N)
        
        # Update rule with proportional threshold
        if k_in[j] > 0:
            new_state = 1 if influence[j] >= w * k_in[j] else 0
            if new_state != state[j]:
                # Recalculate full influence for simplicity and correctness
                state[j] = new_state
                influence = np.zeros(N, dtype=np.int32)
                for i in range(len(src)):
                    if state[src[i]] == 1:
                        influence[dst[i]] += wts[i]
                        
    return evol_active, evol_on_vuln, evol_off_vuln

# =============================================================================
# MAIN SIMULATION WORKFLOW
# =============================================================================

def run_analysis_for_ic(p_initial, N, k_in, w, src, dst, wts):
    """
    Runs all realizations for a single initial condition percentage.
    """
    print(f"  -> Running {REALIZATIONS} realizations for p_initial = {p_initial:.2f}...")
    
    total_steps = NUM_EPOCHS * N
    
    # Arrays to store results from all realizations
    active_trajectories = np.empty((REALIZATIONS, total_steps))
    on_vuln_trajectories = np.empty((REALIZATIONS, total_steps))
    off_vuln_trajectories = np.empty((REALIZATIONS, total_steps))
    
    for r in tqdm(range(REALIZATIONS), desc=f"p_initial={p_initial:.2f}", leave=False):
        # Create a new random initial state for each realization
        init_state = np.random.choice([0, 1], size=N, p=[1 - p_initial, p_initial]).astype(np.int8)
        
        # Run simulation and get the trajectories of metrics
        active_traj, on_vuln_traj, off_vuln_traj = run_proportional_threshold_simulation_jit(
            init_state, k_in, w, total_steps, src, dst, wts
        )
        
        active_trajectories[r, :] = active_traj
        on_vuln_trajectories[r, :] = on_vuln_traj
        off_vuln_trajectories[r, :] = off_vuln_traj
        
    # Calculate mean and std dev across all realizations
    return {
        'm_act': np.mean(active_trajectories, axis=0),
        's_act': np.std(active_trajectories, axis=0),
        'm_on': np.mean(on_vuln_trajectories, axis=0),
        's_on': np.std(on_vuln_trajectories, axis=0),
        'm_off': np.mean(off_vuln_trajectories, axis=0),
        's_off': np.std(off_vuln_trajectories, axis=0)
    }

# =============================================================================
# SCRIPT EXECUTION
# =============================================================================

if __name__ == "__main__":
    start_time = time.time()
    
    # 1. Load and prepare all data structures
    try:
        N_nodes, k_in_arr, src_arr, dst_arr, wts_arr = load_and_prepare_data()
    except FileNotFoundError as e:
        print(f"Fatal Error: {e}")
        exit()

    print(f"\n--- 2. Starting simulation with proportional threshold (w = {PROPORTIONALITY_CONSTANT_W}) ---")
    
    # 2. Iterate through each initial condition, run simulations, and save results
    for p_ic in INITIAL_ACTIVE_FRACTIONS:
        results = run_analysis_for_ic(p_ic, N_nodes, k_in_arr, PROPORTIONALITY_CONSTANT_W, src_arr, dst_arr, wts_arr)
        
        # Save the results for this specific p_ic to a .npz file
        output_filename = OUTPUT_DIR / f"proportional_threshold_w{PROPORTIONALITY_CONSTANT_W}_p{p_ic:.1f}.npz"
        np.savez(output_filename, **results)
        print(f"  -> Results saved to '{output_filename}'")
        
    print(f"\nTotal execution time: {(time.time() - start_time) / 60:.2f} minutes.")