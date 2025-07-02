'''
Author: Jaime Martínez Cazón

Description:
This script performs a verification analysis of the final attractors of the
Gene Regulatory Network under a specific environmental condition (e.g., starvation).
It simulates a large number of individual network realizations starting from
random initial states. Each realization is evolved until it reaches a stable
fixed-point attractor, which is robustly verified using a synchronous update
test.

The primary goal is to determine if all simulations converge to a single
attractor or if multiple stable states exist. The script reports the
distribution of final states, identified by their Hamming distance to a
pre-computed reference attractor.
'''

import time
import numpy as np
import pandas as pd
from numba import njit, prange
from pathlib import Path
from collections import Counter
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
ATTRACTOR_DATA_DIR = BASE_DATA_PATH / "Attractors" # Assumes attractors are in a subfolder

# =============================================================================
# DATA LOADING
# =============================================================================

def load_data_and_conditions():
    """Loads all necessary data and defines simulation conditions."""
    print("--- 1. Loading network, TFs, and reference attractors ---")
    
    # Load network data
    edge_list = pd.read_csv(BASE_DATA_PATH / "giantC_edge_list.csv")
    nodes_id = pd.read_csv(BASE_DATA_PATH / "giantC_nodes_id.csv")
    
    id_to_index = {node_id: idx for idx, node_id in enumerate(nodes_id['ID'].values)}
    N = len(id_to_index)
    
    src = np.array([id_to_index[i] for i in edge_list['Node 1']], dtype=np.int32)
    dst = np.array([id_to_index[j] for j in edge_list['Node 2']], dtype=np.int32)
    wts = edge_list['Weight'].astype(np.int8).to_numpy()

    k_in = np.zeros(N, dtype=np.int32)
    for target in dst:
        k_in[target] += 1
        
    non_fixed_indices = np.where(k_in > 0)[0]
    
    # Define Starvation TFs and load reference attractors
    name_to_id = dict(zip(nodes_id['Gene'], nodes_id['ID']))
    starvation_tfs = ["YMR037C", "YGL073W", "YOL116W", "YNL027W", "YOR028C"] # Add full list
    active_indices = np.array([id_to_index[name_to_id[name]] for name in starvation_tfs if name in name_to_id])
    
    ref_attractors = {}
    for t in THRESHOLDS_TO_RUN:
        path = ATTRACTOR_DATA_DIR / f"reference_attractor_Starvation_t{t}.npz"
        if not path.exists():
            raise FileNotFoundError(f"Reference attractor file not found: {path}")
        ref_attractors[t] = np.load(path)['attractor_state']
        
    print("--- Data loading complete ---")
    return N, k_in, src, dst, wts, non_fixed_indices, active_indices, ref_attractors

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
    
    for attempt in range(max_attempts + 1):
        steps_to_run = base_steps if attempt == 0 else extra_steps
        
        # Asynchronous simulation for a block of steps
        influence = np.zeros(N, dtype=np.int32)
        for i in range(len(src)):
            if state[src[i]] == 1:
                influence[dst[i]] += wts[i]
                
        for _ in range(steps_to_run):
            j = np.random.randint(N)
            if k_in[j] > 0:
                new_state = 1 if influence[j] >= threshold else 0
                if new_state != state[j]:
                    # This simplified influence update is faster if k_out is small
                    # A full re-calculation or outgoing CSR would be more robust
                    state[j] = new_state
                    # NOTE: A full influence update is required for correctness.
                    # Recalculating here for simplicity.
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
            return state, True  # Return stable state and True
            
    return state, False # Return last state and False if not stabilized

# =============================================================================
# MAIN VERIFICATION WORKFLOW
# =============================================================================

def run_verification_simulation(N, k_in, src, dst, wts, non_fixed_indices,
                                active_indices, ref_attractor, threshold):
    """
    Runs the verification for a specific threshold, simulating all realizations.
    """
    print(f"    -> Generating initial population for {REALIZATIONS} realizations...")
    initial_population = np.zeros((REALIZATIONS, N), dtype=np.int8)
    for r in range(REALIZATIONS):
        nodes_to_activate = np.random.choice(non_fixed_indices, 
                                             size=int(len(non_fixed_indices) * INITIAL_ACTIVE_FRACTION), 
                                             replace=False)
        initial_population[r, nodes_to_activate] = 1
        initial_population[r, active_indices] = 1

    final_distances = np.zeros(REALIZATIONS, dtype=np.int32)
    unstable_count = 0
    steps_per_epoch = N
    
    print(f"    -> Running simulations...")
    for r in tqdm(range(REALIZATIONS), desc=f"Verifying t={threshold}"):
        final_state, is_stable = find_final_attractor_jit(
            initial_population[r], k_in, threshold, src, dst, wts,
            NUM_EPOCHS_BASE * steps_per_epoch,
            EXTRA_EPOCHS_PER_ATTEMPT * steps_per_epoch,
            MAX_CONVERGENCE_ATTEMPTS
        )
        if not is_stable:
            unstable_count += 1
        
        final_distances[r] = np.sum(final_state != ref_attractor)
        
    return final_distances, unstable_count

# =============================================================================
# SCRIPT EXECUTION
# =============================================================================

if __name__ == "__main__":
    total_start_time = time.time()
    
    # 1. Load all data
    try:
        N, k_in, src, dst, wts, non_fixed, active_idx, ref_attractors_map = load_data_and_conditions()
    except FileNotFoundError as e:
        print(f"Fatal Error: {e}")
        exit()

    print("\n--- 2. Starting Attractor Verification for Starvation Condition ---")
    
    # 2. Run verification for each threshold
    for t_val in THRESHOLDS_TO_RUN:
        print("-" * 60)
        print(f"--- VERIFYING FOR THRESHOLD t={t_val} ---")
        
        ref_att = ref_attractors_map[t_val]
        distances, num_unstable = run_verification_simulation(
            N, k_in, src, dst, wts, non_fixed, active_idx, ref_att, t_val
        )
        
        # 3. Report the results
        print("\n  --- Final Attractor Distribution Report ---")
        attractor_counts = Counter(distances)
        if len(attractor_counts) == 1:
            dist = list(attractor_counts.keys())[0]
            print(f"  Result: All {REALIZATIONS} realizations converged to a single attractor type.")
            print(f"  Hamming distance to reference attractor: {dist}")
        else:
            print(f"  Result: Multiple final attractor types detected ({len(attractor_counts)} distinct distances).")
            for dist, count in sorted(attractor_counts.items()):
                percentage = (count / REALIZATIONS) * 100
                print(f"    - Hamming distance {dist}: {count} realizations ({percentage:.4f}%)")
        
        if num_unstable > 0:
            print(f"  [WARNING] {num_unstable}/{REALIZATIONS} realizations did not converge to a stable fixed point.")

    print(f"\nTotal execution time: {(time.time() - total_start_time) / 60:.2f} minutes.")