'''
Author: Jaime Martínez Cazón

Description:
This script visualizes the convergence trajectories of the network dynamics by
plotting the distribution of Hamming distances over time. It loads pre-computed
simulation data and generates "split violin plots" for different environmental
conditions (YPD vs. Starvation).

For each specified epoch, it displays the probability density of the Hamming
distances. If the variance of the distribution is negligible, it plots a single
point at the median instead, providing a clear and informative summary of the
system's convergence behavior.
'''

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.stats import gaussian_kde
from matplotlib.patches import Patch

# =============================================================================
# SETUP: PLOTTING PARAMETERS AND PATHS
# =============================================================================
THRESHOLDS_TO_PLOT = [2, 3]
EPOCHS_TO_SHOW = [0, 1, 2, 3, 4, 5, 6, 8, 10, 20]
COLOR_YPD = '#17becf'
COLOR_STARVATION = '#bcbd22'
VIOLIN_WIDTH = 0.8

# Setup paths relative to the script's location
try:
    script_dir = Path(__file__).resolve().parent
except NameError:
    script_dir = Path('.').resolve() # Fallback for interactive environments

DATA_DIR = script_dir / "violin_plot_data"
PLOT_OUTPUT_DIR = script_dir / "plots_final_point"
PLOT_OUTPUT_DIR.mkdir(exist_ok=True)

# =============================================================================
# DATA LOADING AND PLOTTING FUNCTIONS
# =============================================================================

def load_data(data_path, threshold):
    """Loads the convergence trajectory data for a specific threshold."""
    filepath = data_path / f"convergence_trajectories_t{threshold}.npz"
    if not filepath.exists():
        print(f"ERROR: Data file not found at '{filepath}'.")
        return None, None
    
    data = np.load(filepath)
    return data['ypd_trajectories'], data['starvation_trajectories']

def plot_half_violin(ax, data_vector, position_index, color, side='left'):
    """
    Plots one half of a violin plot using Kernel Density Estimation (KDE).
    If variance is too low, plots a single point at the median instead.
    """
    if np.std(data_vector) < 1.0:
        # If variance is negligible, just plot the median point
        median_val = np.median(data_vector)
        ax.scatter(position_index, median_val, color=color, s=200, zorder=3, edgecolors='black')
        return

    try:
        kde = gaussian_kde(data_vector)
        min_val, max_val = data_vector.min(), data_vector.max()
        y_axis = np.linspace(min_val, max_val, 200)
        density = kde(y_axis)
        scaled_density = (density / density.max()) * (VIOLIN_WIDTH / 2)
        
        if side == 'left':
            ax.plot(position_index - scaled_density, y_axis, color=color, lw=1.5)
            ax.fill_betweenx(y_axis, position_index - scaled_density, position_index, color=color, alpha=0.5)
        else: # side == 'right'
            ax.plot(position_index + scaled_density, y_axis, color=color, lw=1.5)
            ax.fill_betweenx(y_axis, position_index, position_index + scaled_density, color=color, alpha=0.5)
    except Exception as e:
        print(f"    - KDE failed, plotting median point instead. Error: {e}")
        ax.scatter(position_index, np.median(data_vector), color=color, s=100, marker='x')

def create_convergence_plot(threshold, ypd_traj, starv_traj):
    """
    Generates and saves the full convergence plot for a given threshold.
    """
    print(f"  Generating plot for Threshold = {threshold}...")
    fig, ax = plt.subplots(figsize=(20, 12))
    
    max_epoch_available = ypd_traj.shape[0] - 1
    valid_epochs = [e for e in EPOCHS_TO_SHOW if e <= max_epoch_available]

    for i, epoch in enumerate(valid_epochs):
        # Plot YPD data on the left half of the violin
        plot_half_violin(ax, ypd_traj[epoch, :], i, COLOR_YPD, side='left')
        
        # Plot Starvation data on the right half of the violin
        plot_half_violin(ax, starv_traj[epoch, :], i, COLOR_STARVATION, side='right')
        
    # --- Aesthetics and Formatting ---
    ax.set_title(f'Network Convergence to Attractor (Threshold = {threshold})', fontsize=30, fontweight='bold', pad=20)
    ax.set_xlabel('Time (Epochs)', fontsize=24, fontweight='bold', labelpad=15)
    ax.set_ylabel('Hamming Distance Distribution', fontsize=24, fontweight='bold', labelpad=15)
    ax.tick_params(axis='both', which='major', labelsize=20)
    ax.set_xticks(range(len(valid_epochs)))
    ax.set_xticklabels(valid_epochs)
    ax.grid(axis='y', linestyle='--', alpha=0.7)
    
    # --- Legend ---
    legend_elements = [Patch(facecolor=COLOR_YPD, label='YPD'),
                       Patch(facecolor=COLOR_STARVATION, label='Dual Starvation')]
    legend = ax.legend(handles=legend_elements, loc='upper right', fontsize=22,
                       title='Environmental Condition', title_fontsize=22, frameon=True)
    plt.setp(legend.get_title(), fontweight='bold')
    
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    
    # --- Save Figure ---
    output_filename = PLOT_OUTPUT_DIR / f"convergence_plot_t{threshold}.pdf"
    plt.savefig(output_filename, dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f"  Plot saved to: {output_filename}")

# =============================================================================
# MAIN EXECUTION
# =============================================================================

if __name__ == '__main__':
    print("--- Starting Plot Generation ---")
    for t in THRESHOLDS_TO_PLOT:
        ypd_data, starv_data = load_data(DATA_DIR, t)
        if ypd_data is not None and starv_data is not None:
            create_convergence_plot(t, ypd_data, starv_data)
    
    print("\nAll plots generated successfully.")