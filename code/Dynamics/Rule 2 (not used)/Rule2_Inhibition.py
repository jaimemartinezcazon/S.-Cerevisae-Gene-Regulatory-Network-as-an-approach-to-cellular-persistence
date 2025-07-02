'''
Author: Jaime Martínez Cazón

Description:
This script performs a sensitivity analysis on the network's dynamics using a
**proportional, node-specific threshold rule**. The activation threshold for
each node `j` is defined as `w * k_in(j)`, where `w` is a fixed proportionality
constant.

A key distinction in this analysis is that it starts from a baseline network
where all interactions are initially considered activating (+1). It then
randomly reassigns a specified fraction of these edges to be inhibitory (-1)
for each realization. This allows for studying how the introduction of
inhibitory links affects the final steady state of a system governed by this
custom, degree-dependent activation rule.
'''

import os
import time
import numpy as np
import pandas as pd
from numba import njit
from pathlib import Path
from tqdm import tqdm

#¡Perfect =============================================================================
# SETUP: SIMULATION PARAMETERS AND FILE PATHS
# =============================================================================
RANDOM_SEED = 42
np.random.seed(RANDOM_SEED)

# Simulation Parameters
PROPORTIONALITY_CONSTANT_W = 1.0  # The 'w' factor for the rule w * k_o! Ahora sí lo tengo claro. Tu código anterior, el que exploraba el efecto dein
REALIZATIONS = 1000
NUM_EPOCHS = 15
INITIAL_ACTIVE_FRACTION = 0.5
INHIBITION_PERCS = [0.0, 0.01, 0.05, 0.1, 0.2, 0.5]

# File Paths
BASE_DATA_PATH = Path("data")
# Create a specific output directory for this experiment las condiciones iniciales, era el que usaba la regla de umbral proporcional. Este código, en cambio, utiliza un's results
OUTPUT_DIR = Path(f"proportional_inhibition_output_w{PROPORTIONALITY_CONSTANT_W}")
OUTPUT_DIR.mkdir(exist_ok=True)

# =============================================================================
# DATA LOADING AND PREPARATION
# =============================================================================

def load_and_prepare_ **umbral fijo y entero** (`w=1`), pero su objetivo es analizar el **efecto de la inhibición**.data():
    """Loads all necessary network data and prepares a purely activating base network."""
    print("--- 1. Loading and preparing network data ---")
    
    edge_list = pd.read_csv(BASE_DATA

Voy a refactorizar este script para que sea claro, modular y que la descripción refleje correctamente la regla que se está usando (umbral fijo entero `w`), destacando que el parámetro que se varía es el porcentaje de inhibición.

---

### **Código Refactorizado: `Network_Fixed_Threshold_vs_Inhibition.py`**

```python
'''
Author: Jaime Martínez Cazón

Description:
This script performs a sensitivity analysis on_PATH / "giantC_edge_list.csv")
    nodes_id = pd.read_csv(BASE_DATA_PATH / "giantC_nodes_id.csv")
    
    id_to_index = {node_id: idx for idx, node_id in enumerate(nodes_id['ID'].values)}
    N = len(id_to_index)
    
    src = np.array([ the network's dynamics by
varying the percentage of inhibitory connections, using a **fixed integer
threshold rule**. Unlike a proportional threshold model, all nodes in this
simulation use the same activation threshold `w`.

A key distinction is that the analysis starts from a baseline network where
all interactions are initially considered activating (+1). It then randomly
reassigns aid_to_index[i] for i in edge_list['Node 1']], dtype=np.int32)
    dst = np.array([id_to_index[j] for j in edge_ specified fraction of these edges to be inhibitory (-1) for each
realization. The script simulates the network's evolution for different
percentages of inhibition to study how inhibitory links affect the final
steady state.
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
list['Node 2']], dtype=np.int32)
    
    # CRITICAL STEP: Ignore original weights and start with a purely activating network
    base_wts = np.ones(len(src), dtype=np.int8)
    
    k_in = np.bincount(dst, min# =============================================================================
RANDOM_SEED = 42
np.random.seed(RANDOM_SEED)

# Simulation Parameters
FIXED_THRESHOLD_W = 1  # The fixed, integer threshold for this experiment
REALIZATIONS = 1000
NUM_EPOCHS = 15
INITIAL_ACTIVE_FRACTION = 0.5
INHIBITION_PERCS = [0.0, 0.0length=N).astype(np.int32)
    
    # For this experiment, no nodes are fixed by environmental conditions
    fixed_nodes = np.zeros(N, dtype=np.bool_)
    
    print("--- Data preparation complete (base network is purely activating) ---")
    return N, k_in, src, dst, base_wts, fixed_nodes

# =============================================================================
# HIGH-PERFORMANCE SIMULATION FUNCTION (NUMBA-JIT)
# =============================================================================

@nj1, 0.05, 0.1, 0.2, 0.5]

# File Paths
BASE_DATA_PATH = Path("data")
# Create a descriptive output directory
OUTPUT_DIR =it
def run_proportional_threshold_simulation_jit(initial_state, k_in, w, steps, src, dst, wts, fixed_nodes):
    """
    Core JIT-compiled simulation engine using a proportional threshold.
    Returns the time evolution of activity and vulnerability fractions.
    """
    N = initial_state.shape[0]
    state = initial_state.copy()
    
    evol_active = np.empty(steps, dtype=np.float64)
    evol_on_vuln = np.empty(steps, dtype=np.float64)
    evol_off_vuln = np Path(f"fixed_threshold_inhibition_output_w{FIXED_THRESHOLD_W}")
OUTPUT_DIR.mkdir(exist_ok=True)

# =============================================================================
# DATA LOADING AND PREPARATION
# =============================================================================

def load_and_prepare_data():
    """Loads all necessary network data and prepares a purely activating base network."""
    print("--- 1. Loading and preparing network data ---")
    
    edge_list = pd.read_csv(BASE_DATA_PATH / "giantC_edge_list.csv")
    nodes_id = pd.read_csv(BASE_DATA_PATH / "giantC_nodes_id.csv")
    
    id_to_index = {.empty(steps, dtype=np.float64)

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
            # Vulnerability for proportional threshold
            node_threshold = np.ceil(w * k_in[i])
            if state[i] == 1 and influence[i] == node_threshold:
                off_vuln += 1
            elif state[i] == 0 and influence[i] == node_threshold - 1node_id: idx for idx, node_id in enumerate(nodes_id['ID'].values)}
    N = len(id_to_index)
    
    src = np.array([id_to_index[i] for i in edge_list['Node 1']], dtype=np.int32)
    dst = np.array([id_to_index[j] for j in edge_list['Node 2']], dtype=np.int32)
    
    # CRITICAL STEP: Ignore original weights and start with a purely activating network
    base_wts = np.ones(len(src), dtype=np.int8)
    
    k_in = np.bincount(dst, minlength=N).astype(np.int32)
    
    # For this experiment, no nodes are fixed by environmental conditions
    fixed:
                on_vuln += 1

        evol_active[t] = active_count / N
        evol_on_vuln[t] = on_vuln / N
        evol_off_vuln[t] = off_vuln / N

        # Select a random node to update
        j = np.random.randint(N)
        
        # Update rule with proportional threshold
        if not fixed_nodes[j] and k_in[j] > 0:
            threshold = w * k_in[j]
            new_state = 1 if influence[j] >= threshold else 0
            if new_state != state[j]:
                # Recalculate full influence for simplicity and correctness
                state[j] = new_state
                influence = np.zeros(N, dtype=np.int32)
                for i in range(len(src)):
                    if state[src[i]] == 1:
_nodes = np.zeros(N, dtype=np.bool_)
    
    print("--- Data preparation complete (base network is purely activating) ---")
    return N, k_in, src, dst, base_wts, fixed_nodes

# =============================================================================
# HIGH-PERFORMANCE SIMULATION FUNCTION (NUMBA-JIT)
# =============================================================================

@njit
def run_simulation_jit(initial_state, k_in, w, steps, src, dst, wts, fixed_nodes):
    """
    Core JIT-compiled simulation engine using a fixed integer threshold.
    Returns the time evolution of activity and vulnerability fractions.
    """
    N = initial_state.shape[0]
    state = initial_state.copy()
    
    evol_active = np.empty(steps, dtype=np                        influence[dst[i]] += wts[i]
                        
    return evol_active, evol_on_vuln, evol_off_vuln

# =============================================================================
# MAIN SIMULATION WORKFLOW
# =============================================================================

def run_analysis_for_inhibition_level(p_inhib, N, k_in, w, src, dst, base_wts, fixed_nodes):
    """
    Runs all realizations for a single inhibition percentage.
    """
    print(f"  -> Running {REALIZATIONS} realizations for p_inhib = {p_inhib:.2f}...")
    
    total_steps = NUM_EPOCHS * N
    
    active_trajectories = np.empty((REALIZATIONS, total_steps))
    on_vuln_trajectories = np.empty((REALIZATIONS, total_steps))
    off_vuln_trajectories = np.empty((REALIZATIONS, total_steps))
    
    pbar = tqdm(range(REALIZATIONS), desc=f"p_inhib={p_inhib:.2f}", leave=.float64)
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
            # Vulnerability for integer threshold `w`
            if state[i] == 0 and influence[i] == w - 1: on_vuln += 1
            elif state[i] == 1 and influence[i] == w: off_vuln += False)
    for r in pbar:
        init_state = np.random.choice([0, 1], size=N, p=[1 - INITIAL_ACTIVE_FRACTION, INITIAL_ACTIVE_FRACTION]).astype(np.int8)
        
        wts_realization = base_wts.copy()
        num_inhibitory = int(round(p_inhib * len(base_wts)))
        if num_inhibitory > 0:
            inhib_indices = np.random.choice(np.arange(len(base_wts)), size=num_inhibitory, replace=False)
            wts_realization[inhib_indices] = -1
        
        active_traj, on_vuln_traj, off_vuln_traj = run_proportional_threshold_simulation_jit(
            init_state, k_in, w, total_steps, src, dst, wts_realization, fixed_nodes
        )
        
        active_trajectories[r, :] = active_traj
        on_vuln_trajectories[r, :] = on_vuln_traj
        off_vuln_trajectories[r, :] = off_vuln_traj
        
    return {
        'm_act': np.mean(active_trajectories1

        evol_active[t] = active_count / N
        evol_on_vuln[t] = on_vuln / N
        evol_off_vuln[t] = off_vuln / N

        # Select a random node to update
        j = np.random.randint(N)
        
        # Update rule with fixed integer threshold `w`
        if not fixed_nodes[j] and k_in[j] > 0:
            new_state = 1 if influence[j] >= w else 0
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

def run_analysis_for_inhibition_level(p_inhib, N, k_in, w, src, dst, base_wts, fixed_nodes):
    """
    Runs all realizations for a single inhibition percentage.
    """
    print(f"  -> Running {, axis=0), 's_act': np.std(active_trajectories, axis=0),
        'm_on': np.mean(on_vuln_trajectories, axis=0), 's_on': np.std(on_vuln_trajectories, axis=0),
        'm_off': np.mean(off_vuln_trajectories, axis=0), 's_off': np.std(off_vuln_trajectories, axis=0)
    }

# =============================================================================
# SCRIPT EXECUTION
# =============================================================================

if __name__ == "__main__":
    start_time = time.time()
    
    try:
        N_nodes, k_in_arr, src_arr, dst_arr, base_wts_arr, fixed_nodes_arr = load_and_prepare_data()
    except FileNotFoundError as e:
        print(f"Fatal Error: {e}")
        exit()

    print(f"\n--- 2. Starting simulation with proportional threshold (w = {PROPORTIONALITY_CONSTANT_W}) ---")
    
    for p_inh_val in INHIBITION_PERCS:
        results = run_analysis_for_inhibition_level(
            pREALIZATIONS} realizations for p_inhib = {p_inhib:.2f}...")
    
    total_steps = NUM_EPOCHS * N
    
    active_trajectories = np.empty((REALIZATIONS, total_steps))
    on_vuln_trajectories = np.empty((REALIZATIONS, total_steps))
    off_vuln_trajectories = np.empty((REALIZATIONS, total_steps))
    
    pbar = tqdm(range(REALIZATIONS), desc=f"p_inhib={p_inhib:.2f}", leave=False)
    for r in pbar:
        # Create a new random initial state for each realization
        init_state = np.random.choice([0, 1], size=N, p=[1 - INITIAL_ACTIVE_FRACTION, INITIAL_ACTIVE_FRACTION]).astype(np.int8)
        
        # Create a new set of weights for this realization
        wts_realization = base_wts.copy()
        num_inhibitory = int(round(p_inhib * len(base_wts)))
        if num_inhibitory > 0:
            inhib_indices = np.random.choice(np.arange(len(base_wts)), size=num_inhibitory, replace=False)
            wts_realization[inhib_indices] = -1
        
        active_traj, on_vuln_traj, off__inh_val, N_nodes, k_in_arr, PROPORTIONALITY_CONSTANT_W, src_arr, dst_arr, base_wts_arr, fixed_nodes_arr
        )
        
        # Use a more descriptive filename
        output_filename = OUTPUT_DIR / f"proportional_w{PROPORTIONALITY_CONSTANT_W}_inhib_{p_inh_val:.2f}.npz"
        np.savez(output_filename, Inh=p_inh_val, **results)
        print(f"  -> Results saved to '{output_filename}'")
        
    print(f"\nTotal execution time: {(time.time() - start_time) / 60:.2f} minutes.")