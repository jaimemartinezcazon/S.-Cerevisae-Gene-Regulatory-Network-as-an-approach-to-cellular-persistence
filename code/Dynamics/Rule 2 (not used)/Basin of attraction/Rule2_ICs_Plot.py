'''
Author: Jaime Martínez Cazón

Description:
This script visualizes the results from the network dynamics simulation that used
a **proportional threshold** (where threshold is proportional to a node's in-degree).
It loads a series of .npz files, each corresponding to a different initial
fraction of active nodes, for a specific proportionality constant `w`.

The script generates a multi-panel plot showing the time evolution of:
- The average fraction of active nodes.
- The average fraction of OFF-vulnerable nodes.
- The average fraction of ON-vulnerable nodes.

Each curve represents a different initial condition, allowing for a clear
visualization of how the starting state affects the network's trajectory under
this specific update rule.
'''

import glob
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

# =============================================================================
# SETUP: PLOTTING PARAMETERS AND PATHS
# =============================================================================
# The proportionality constant 'w' for which to generate plots
PROPORTIONALITY_CONSTANT_W = 0.2

# Path to the directory containing the simulation results
DATA_DIR = Path("proportional_threshold_output")

# =============================================================================
# DATA LOADING AND PLOTTING FUNCTIONS
# =============================================================================

def load_simulation_data_for_w(data_dir, w):
    """
    Loads all .npz files for a given proportionality constant 'w'.

    Args:
        data_dir (Path): The directory containing the data files.
        w (float): The proportionality constant to load data for.

    Returns:
        list: A list of dictionaries, where each dictionary contains the data
              for one initial condition. Returns an empty list if no files are found.
    """
    print(f"--- 1. Loading data for w = {w} from '{data_dir}' ---")
    
    # Find all files matching the pattern for the given w
    file_pattern = f"proportional_threshold_w{w}_p*.npz"
    file_paths = sorted(data_dir.glob(file_pattern))
    
    if not file_paths:
        print(f"Error: No data files found for w={w} in the specified directory.")
        return []

    loaded_data = []
    for fpath in file_paths:
        try:
            data = np.load(fpath)
            # Store all relevant arrays and parameters in a dictionary
            loaded_data.append({
                'p_initial': float(data['IC']),
                'm_act': data['m_act'], 's_act': data['s_act'],
                'm_off': data['m_off'], 's_off': data['s_off'],
                'm_on': data['m_on'], 's_on': data['s_on']
            })
        except Exception as e:
            print(f"Warning: Could not load or parse file '{fpath.name}'. Skipping. Error: {e}")
            
    return loaded_data


def create_trajectory_plot(simulation_data, w):
    """
    Generates the 3-panel plot showing the time evolution of network metrics.
    """
    if not simulation_data:
        print("No data available to generate plot.")
        return

    print("--- 2. Generating trajectory plot ---")
    
    fig, axes = plt.subplots(3, 1, figsize=(12, 18), sharex=True)
    fig.suptitle(f'Network Dynamics with Proportional Threshold (w = {w})', 
                 fontsize=28, fontweight='bold', y=0.98)
    
    # Get the number of steps from the first loaded data entry
    num_steps = len(simulation_data[0]['m_act'])
    x_axis = np.arange(num_steps)

    plot_titles = ['Active Node Fraction', 'OFF-Vulnerable Fraction', 'ON-Vulnerable Fraction']
    metric_keys = [('m_act', 's_act'), ('m_off', 's_off'), ('m_on', 's_on')]

    for ax, title, (mean_key, std_key) in zip(axes, plot_titles, metric_keys):
        ax.set_title(title, fontsize=24, fontweight='bold')
        
        for data_point in simulation_data:
            p_initial = data_point['p_initial']
            mean_traj = data_point[mean_key]
            std_traj = data_point[std_key]
            
            line, = ax.plot(x_axis, mean_traj, lw=2, label=f'{p_initial*100:.0f}%')
            ax.fill_between(x_axis, mean_traj - std_traj, mean_traj + std_traj,
                            color=line.get_color(), alpha=0.2)
        
        # Formatting for each subplot
        ax.set_ylim(0, 1.05)
        ax.set_xlim(0, num_steps)
        ax.tick_params(labelsize=18)
        ax.grid(True, linestyle=':', alpha=0.6)

    # Set shared labels
    axes[1].set_ylabel('Fraction', fontsize=24, fontweight='bold')
    axes[-1].set_xlabel('Time (steps)', fontsize=24, fontweight='bold')

    # Create a single legend for the entire figure
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, title='Initial Active %', loc='center right',
               fontsize=20, title_fontsize=22)
    
    plt.tight_layout(rect=[0, 0, 0.85, 0.95]) # Adjust rect to make space for legend
    plt.show()

# =============================================================================
# MAIN EXECUTION
# =============================================================================

if __name__ == "__main__":
    # 1. Load all simulation data for the specified 'w'
    all_data = load_simulation_data_for_w(DATA_DIR, PROPORTIONALITY_CONSTANT_W)
    
    # 2. Generate the plot if data was successfully loaded
    if all_data:
        create_trajectory_plot(all_data, PROPORTIONALITY_CONSTANT_W)
    else:
        print("Plot generation skipped due to data loading errors.")