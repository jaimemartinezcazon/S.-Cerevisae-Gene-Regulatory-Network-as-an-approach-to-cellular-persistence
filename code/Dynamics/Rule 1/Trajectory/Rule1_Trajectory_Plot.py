'''
Author: Jaime Martínez Cazón

Description:
This script visualizes the results from the network dynamics simulation. It
loads the pre-computed trajectory data from a .npz file and generates a
multi-panel plot showing the time evolution of several key metrics:
- Overall network activity (fraction of active nodes)
- Candidate gene set activity
- On-vulnerable and Off-vulnerable node fractions

The plot is designed with a vertical layout, a shared x-axis in units of
epochs, and specific styling to ensure clarity and academic presentation quality.
'''

import os
import numpy as np
import matplotlib.pyplot as plt

# =============================================================================
# SETUP: FILE PATHS AND PLOT CONFIGURATION
# =============================================================================
INPUT_FILENAME = 'simulation_trajectory_data.npz'
END_EPOCH_PLOT = 45  # The upper limit for the x-axis in epochs

# =============================================================================
# DATA LOADING FUNCTION
# =============================================================================

def load_simulation_data(filepath):
    """
    Loads the simulation results from a .npz file.

    Args:
        filepath (str): The path to the input .npz file.

    Returns:
        dict or None: A dictionary containing the loaded data, or None if an error occurs.
    """
    if not os.path.exists(filepath):
        print(f"Error: Simulation data file not found: '{filepath}'")
        print("Please run the network dynamics simulation script first.")
        return None
    
    print(f"Loading simulation data from '{filepath}'...")
    try:
        data = np.load(filepath)
        # Check for essential keys to ensure the file is valid
        if 'thresholds' not in data or 'total_steps' not in data:
            raise KeyError("Essential keys ('thresholds', 'total_steps') are missing from the data file.")
        return data
    except Exception as e:
        print(f"Error loading data from '{filepath}': {e}")
        return None

# =============================================================================
# PLOTTING FUNCTIONS
# =============================================================================

def plot_metric_trajectory(ax, data, metric_key, title, colors):
    """
    Plots the time evolution for a single metric (e.g., activity) for all thresholds.
    """
    thresholds = data['thresholds']
    total_steps = int(data['total_steps'])
    x_plot = np.arange(total_steps)
    
    for w in thresholds:
        mean_key = f'mean_{metric_key}_{w}'
        std_key = f'std_{metric_key}_{w}'
        
        if mean_key in data and std_key in data:
            mean_values = data[mean_key]
            std_values = data[std_key]
            color = colors.get(w, 'gray')
            
            ax.plot(x_plot, mean_values, lw=2, color=color, label=f't={w}')
            ax.fill_between(x_plot, mean_values - std_values, mean_values + std_values,
                            color=color, alpha=0.3)
    ax.set_ylim(0, 1.05)
    ax.set_ylabel(title, fontsize=24, fontweight='bold')


def configure_plot_aesthetics(fig, axes, data):
    """Applies final aesthetic configurations to the entire figure."""
    # --- Calculate tick positions ---
    steps_per_epoch = int(data['steps_per_epoch'])
    epochs_per_segment = int(data['epochs_per_segment'])
    num_segments = int(data['num_segments'])
    
    end_step_plot = END_EPOCH_PLOT * steps_per_epoch
    transition_epochs = [epochs_per_segment * i for i in range(1, num_segments)]
    tick_steps = [0] + [steps_per_epoch * e for e in transition_epochs] + [end_step_plot]
    tick_labels = [0] + transition_epochs + [END_EPOCH_PLOT]
    
    # --- Configure all axes ---
    for i, ax in enumerate(axes):
        ax.set_xlim(0, end_step_plot)
        ax.set_xticks(tick_steps)
        ax.set_yticks([0, 0.5, 1])
        ax.grid(True, linestyle=':', alpha=0.6)
        ax.tick_params(axis='y', labelsize=30)
        
        # Draw vertical lines for transitions
        for epoch in transition_epochs:
            ax.axvline(epoch * steps_per_epoch, color='gray', linestyle='--', lw=1.5)
        
        # Configure x-axis labels only for the bottom plot
        if i == len(axes) - 1:
            ax.set_xticklabels(tick_labels)
            ax.tick_params(axis='x', labelsize=30)
            ax.set_xlabel('Time (epochs)', fontsize=42, fontweight='bold')
        else:
            ax.tick_params(axis='x', labelbottom=False)
            
    # --- Configure Figure-level elements ---
    fig.supylabel('Active Nodes Ratio', x=0.04, fontsize=42, fontweight='bold')
    
    handles, labels = axes[0].get_legend_handles_labels()
    legend = fig.legend(handles, labels, title='Threshold', loc='upper right', 
                        bbox_to_anchor=(0.98, 0.9), frameon=False, fontsize=36,
                        title_fontsize=36)
    plt.setp(legend.get_title(), fontweight='bold')
    
    plt.tight_layout(rect=[0.05, 0, 0.85, 1]) # Adjust rect to make space for labels/legend

# =============================================================================
# MAIN EXECUTION
# =============================================================================

if __name__ == "__main__":
    # 1. Load the simulation data
    simulation_data = load_simulation_data(INPUT_FILENAME)
    
    if simulation_data:
        # 2. Setup the plot figure
        fig, axes = plt.subplots(4, 1, figsize=(12, 20), sharex=True)
        
        plot_configs = [
            {'key': 'active', 'title': 'Overall Activity', 'ax': axes[0]},
            {'key': 'cand', 'title': 'Candidate Set Activity', 'ax': axes[1]},
            {'key': 'onv', 'title': 'On-Vulnerable Fraction', 'ax': axes[2]},
            {'key': 'offv', 'title': 'Off-Vulnerable Fraction', 'ax': axes[3]}
        ]
        
        threshold_colors = {2: '#bcbd22', 3: '#ff7f0e'}
        
        # 3. Create each subplot
        print("Generating plots...")
        for config in plot_configs:
            plot_metric_trajectory(
                config['ax'],
                simulation_data,
                config['key'],
                config['title'],
                threshold_colors
            )
            
        # 4. Apply final aesthetic configurations
        configure_plot_aesthetics(fig, axes, simulation_data)
        
        print("Displaying final plot.")
        plt.show()