#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# Load the CSV file into a DataFrame
df = pd.read_csv('increasing.csv')

# Extract the relevant columns for plotting
haplotypes = df["Read"].tolist()

# Function to extract values from each tuple string
def extract_edit_distance(value):
    return int(value.split(', ')[-1].strip(')'))

def extract_runtime(value):
    return float(value.split(', ')[0].strip('('))

def extract_memory(value):
    return float(value.split(', ')[1].strip())

# Extract data for each coverage level
edit_distances_13M = df["13M"].apply(extract_edit_distance).tolist()
edit_distances_25M = df["25M"].apply(extract_edit_distance).tolist()
edit_distances_49M = df["49M"].apply(extract_edit_distance).tolist()

runtime_13M = df["13M"].apply(extract_runtime).tolist()
runtime_25M = df["25M"].apply(extract_runtime).tolist()
runtime_49M = df["49M"].apply(extract_runtime).tolist()

memory_13M = df["13M"].apply(extract_memory).tolist()
memory_25M = df["25M"].apply(extract_memory).tolist()
memory_49M = df["49M"].apply(extract_memory).tolist()

# Initialize a figure with 3 horizontally stacked subplots
fig, axes = plt.subplots(1, 3, figsize=(10, 2.5))

# X-axis positions
x = np.arange(len(haplotypes))
width = 0.2  # Width of bars

# Plot 1: Edit Distance (Log-Scaled)
# Apply log transformation to edit distances
edit_distances_13M_log = np.log10(edit_distances_13M)
edit_distances_25M_log = np.log10(edit_distances_25M)
edit_distances_49M_log = np.log10(edit_distances_49M)

axes[0].bar(x - width, edit_distances_13M_log, width, label='13H', zorder=3, color=plt.get_cmap('tab10')(0))
axes[0].bar(x, edit_distances_25M_log, width, label='25H', zorder=3, color=plt.get_cmap('tab10')(1))
axes[0].bar(x + width, edit_distances_49M_log, width, label='49H', zorder=3, color=plt.get_cmap('tab10')(2))

axes[0].set_xlabel('Haplotypes', fontsize=12)
axes[0].set_ylabel('Edit Distance', fontsize=12)
axes[0].set_xticks(x)
axes[0].set_xticklabels(haplotypes, fontsize=11)
y_ticks = np.arange(0, 5)
axes[0].set_yticks(y_ticks)
axes[0].set_yticklabels([f'$10^{int(tick)}$' for tick in y_ticks], fontsize=12)
axes[0].grid(axis='y', linestyle='--', alpha=0.6, zorder=0)

# Plot 2: Runtime (hours)
# Convert runtime from seconds to hours
runtime_13M_hours = np.array(runtime_13M) / 3600
runtime_25M_hours = np.array(runtime_25M) / 3600
runtime_49M_hours = np.array(runtime_49M) / 3600

axes[1].bar(x - width, runtime_13M_hours, width, label='13H', zorder=3, color=plt.get_cmap('tab10')(0))
axes[1].bar(x, runtime_25M_hours, width, label='25H', zorder=3, color=plt.get_cmap('tab10')(1))
axes[1].bar(x + width, runtime_49M_hours, width, label='49H', zorder=3, color=plt.get_cmap('tab10')(2))

axes[1].set_xlabel('Haplotypes', fontsize=12)
axes[1].set_ylabel('Runtime (hours)', fontsize=12)
axes[1].set_xticks(x)
axes[1].set_xticklabels(haplotypes, fontsize=12)
axes[1].grid(axis='y', linestyle='--', alpha=0.6, zorder=0)
axes[1].tick_params(axis='y', labelsize=12)

# Plot 3: Memory Usage (GB)
axes[2].bar(x - width, memory_13M, width, label='13H', zorder=3, color=plt.get_cmap('tab10')(0))
axes[2].bar(x, memory_25M, width, label='25H', zorder=3, color=plt.get_cmap('tab10')(1))
axes[2].bar(x + width, memory_49M, width, label='49H', zorder=3, color=plt.get_cmap('tab10')(2))

axes[2].set_xlabel('Haplotypes', fontsize=12)
axes[2].set_ylabel('Memory Usage (GB)', fontsize=12)
axes[2].set_xticks(x)
axes[2].set_xticklabels(haplotypes, fontsize=12)
axes[2].grid(axis='y', linestyle='--', alpha=0.6, zorder=0)
axes[2].tick_params(axis='y', labelsize=12)

# Adding a legend to the first plot
handles, labels = axes[0].get_legend_handles_labels()
fig.legend(handles, labels, loc='upper center', ncol=3, bbox_to_anchor=(0.5, 1.12), fontsize=12)

# Adjust layout and add space between subplots
plt.tight_layout()
plt.subplots_adjust(wspace=0.4)  # Add space between the plots

# Save the figure
plt.savefig('increasing.pdf', format='pdf', dpi=1200, bbox_inches='tight')