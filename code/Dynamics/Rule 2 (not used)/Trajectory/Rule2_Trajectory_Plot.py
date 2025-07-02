'''
Author: Jaime Martínez Cazón

Description:
This script visualizes the results from the network dynamics simulation that
used a **proportional, node-specific threshold rule** (`threshold = w * k_in(j)`).
It loads a series of .npz files, each corresponding to a different
proportionality constant `w`, while the initial condition was held fixed.

The script generates a multi-panel plot showing the time evolution of network
metrics (activity, vulnerability). Each curve in the plot represents a
different value of `w`, allowing for a clear visualization of how this key
parameter affects the network's trajectory and final steady state.
'''

import glob
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

# =============================================================================
# SETUP: FILE PATHS AND PARAMETERS
# =============================================================================
# Path to the directory containing the simulation results
DATA_DIR = Path("proportional_threshold_vs_w_output")

# =============================================================================
# DATA LOADING AND PLOTTING FUNCTIONS
# =============================================================================

def load_simulation_data_for_w_range(data_dir):
    """
    Loads all .npz files from the proportionality factor experiment.

    Args:
        data_dir (Path): The directory containing the data files.

    Returns:
        list: A list of dictionaries, where each dictionary contains the data
              for one value of 'w'. Returns an empty list if no files found.
    """
    print(f"--- 1. Loading data from '{data_dir}' ---")
    
    file_pattern = "proportional_threshold*.npz"
    file_paths = sorted(data_dir.glob(file_pattern))
    
    if not file_paths:
        print(f"Error: No data files found in the specified directory.")
        return []

    loaded_data = []
    for fpath in file_paths:
        try:
            data = np.load(fpath)
            # The key 'w' was used in the simulation script
            loaded_data.append({
                'w_factor': float(data['w']),
                'm_act': data['m_act'], 's_act': data['s_act'],
                'm_off': data['m_off'], 's_off': data['s_off'],
                'm_on': data['m_on'], 's_on': data['s_on']
            })
        except Exception as e:
            print(f"Warning: Could not load or parse file '{fpath.name}'. Skipping. Error: {e}")
            
    # Sort data by 'w' factor for a clean legend
    loaded_data.sort(key=lambda x: x['w_factor'])
    return loaded_data


def create_trajectory_plot(simulation_data):
    """
    Generates the 3-panel plot showing the time evolution of network metrics.
    """
    if not simulation_data:
        print("No data available to generate plot.")
        return

    print("--- 2. Generating trajectory plot ---")
    
    fig, axes = plt.subplots(3, 1, figsize=(12, 18), sharex=True)
    fig.suptitle('Network Dynamics vs. Proportionality Factor (w)', 
                 fontsize=28, fontweight='bold', y=0.98)
    
    num_steps = len(simulation_data[0]['m_act'])
    x_axis = np.arange(num_steps)

    plot_titles = ['Active Node Fraction', 'OFF-Vulnerable Fraction', 'ON-Vulnerable Fraction']
    metric_keys = [('m_act', 's_act'), ('m_off', 's_off'), ('m_on', 's_on')]

    for ax, title, (mean_key, std_key) in zip(axes, plot_titles, metric_keys):
        ax.set_title(title, fontsize=24, fontweight='bold')
        
        for data_point in simulation_data:
            w_val = data_point['w_factor']
            mean_traj = data_point[mean_key]
            std_traj = data_point[std_key]
            
            line, = ax.plot(x_axis, mean_traj, lw=2, label=f'w={w_val:.1f}')
            ax.fill_between(x_axis, mean_traj - std_traj, mean_traj + std_traj,
                            color=line.get_color(), alpha=0.2)
        
        ax.set_ylim(0, 1.05)
        ax.set_xlim(0, num_steps)
        ax.tick_params(labelsize=18)
        ax.grid(True, linestyle=':', alpha=0.6)

    axes[1].set_ylabel('Fraction', fontsize=24, fontweight='bold')
    axes[-1].set_xlabel('Time (steps)', fontsize=24, fontweight='bold')

    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, title='Factor (w)', loc='center right',
               fontsize=20, title_fontsize=22)
    
    plt.tight_layout(rect=[0, 0, 0.85, 0.95])
    plt.show()

# =============================================================================
# MAIN EXECUTION
# =============================================================================

if __name__ == "__main__":
    # 1. Load all simulation data
    all_data = load_simulation_data_for_w_range(DATA_DIR)
    
    # 2. Generate the plot if data was successfully loaded
    if all_data:
        create_trajectory_plot(all_data)
    else:
        print("Plot generation skipped due to data loading errors.")