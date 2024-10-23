import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Load the CSV files from current directory paths
df_pan_genie = pd.read_csv('PanGenie.csv', index_col=0)
df_phi = pd.read_csv('PHI_MIQP.csv', index_col=0)
df_vg = pd.read_csv('VG.csv', index_col=0)

# Coverage levels and tools
coverages = ['0.1x', '0.5x', '1x', '2x', '5x', '10x', '15x']
last_cov = ['16.26x', '12.91x', '18.20x', '12.85x', '15.04x']

# Function to extract the edit distance for a given tool, coverage, and read
def get_edit_distance(df, read, coverage):
    return int(df.loc[read, coverage].split(', ')[-1][:-1])

# Plotting bar plots
fig, axes = plt.subplots(nrows=5, ncols=1, figsize=(8, 9))

reads = df_pan_genie.index

# Define y-ticks to ensure all powers of 10 are included
y_ticks = np.arange(0, 7, 2)

# Alphabet labels for subplots
alphabet_labels = ['(A)', '(B)', '(C)', '(D)', '(E)']

# Plot each read separately
for i, read in enumerate(reads):
    # Extracting edit distances for each coverage level and tool
    edit_distances_phi = [get_edit_distance(df_phi, read, coverage) for coverage in coverages]
    edit_distances_vg = [get_edit_distance(df_vg, read, coverage) for coverage in coverages]
    edit_distances_pan_genie = [get_edit_distance(df_pan_genie, read, coverage) for coverage in coverages]
    
    # Convert to log10 scale
    edit_distances_phi_log = np.log10(edit_distances_phi)
    edit_distances_vg_log = np.log10(edit_distances_vg)
    edit_distances_pan_genie_log = np.log10(edit_distances_pan_genie)
    
    # X-axis positions
    x = np.arange(len(coverages))
    width = 0.2  # Adjusted bar width for the remaining plots
    
    # Plot bars for each tool at different positions in the required order
    axes[i].bar(x - width, edit_distances_phi_log, width, label='PHI', zorder=3)
    axes[i].bar(x, edit_distances_vg_log, width, label='VG', zorder=3)
    axes[i].bar(x + width, edit_distances_pan_genie_log, width, label='PanGenie', zorder=3)
    
    axes[i].set_ylabel('Edit\ndistance', fontsize=13)
    axes[i].set_xticks(x)
    
    # Adjust the x-tick labels for the last coverage
    new_cov = [last_cov[i] if cov == '15x' else cov for cov in coverages]
    axes[i].set_xticklabels(new_cov, fontsize=13)
    axes[i].set_yticks(y_ticks)
    axes[i].set_yticklabels([f'$10^{int(tick)}$' for tick in y_ticks], fontsize=13)
    axes[i].set_xlabel('Coverage', fontsize=13)
    axes[i].grid(axis='y', linestyle='--', alpha=0.6, zorder=0)
    
    # Add the title inside the bar plot
    # axes[i].set_title(f'{alphabet_labels[i]} {read}', fontsize=12, loc='center')
    axes[i].set_title(f'{alphabet_labels[i]} {read}', fontsize=13, loc='center', fontweight='bold')

axes[-1].set_xlabel('Coverage', fontsize=13)

# Adding the legend at the top of the first plot
handles, labels = axes[0].get_legend_handles_labels()
fig.legend(handles, labels, loc='upper right', fontsize=14, bbox_to_anchor=(1.23, 0.5), title='Method', title_fontsize='14')

# Adjust layout to add space between plots and ensure labels are visible
plt.subplots_adjust(hspace=0.5)
plt.tight_layout(rect=[0, 0, 1, 0.98])

plt.savefig('edit_distances.pdf', format='pdf', bbox_inches='tight', dpi=1200)