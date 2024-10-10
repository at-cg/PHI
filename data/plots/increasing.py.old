#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# Load the CSV file into a DataFrame
df_updated = pd.read_csv('increasing.csv')

# Function to extract values from each tuple string
def extract_edit_distance(value):
    return int(value.split(', ')[-1].strip(')'))

def extract_runtime(value):
    return float(value.split(', ')[0].strip('('))

def extract_memory(value):
    return float(value.split(', ')[1].strip())

# Add columns for 3M, 7M, 13M, 25M, and 49M with extracted data
coverage_levels = ["3H", "7H", "13H", "25H", "49H"]

for cov in coverage_levels:
    df_updated[cov] = df_updated[cov].apply(lambda x: tuple(map(float, x.strip('()').split(', '))))

# Extract the relevant columns for plotting
haplotypes = df_updated["Read"].tolist()

edit_distances = {cov: df_updated[cov].apply(lambda x: x[2]).tolist() for cov in coverage_levels}
runtimes = {cov: df_updated[cov].apply(lambda x: x[0] / 3600).tolist() for cov in coverage_levels}  # Convert to hours
memory_usages = {cov: df_updated[cov].apply(lambda x: x[1]).tolist() for cov in coverage_levels}

x = np.arange(len(haplotypes))
width = 0.15  # Width of bars for edit distance plot

# Plot 1: Edit Distance (Log Scale)
plt.figure(figsize=(8, 3))
for i, cov in enumerate(coverage_levels):
    plt.bar(x - (2-i)*width, np.log10(edit_distances[cov]), width, label=f'{cov}', color=plt.get_cmap('tab10')(i), zorder=3)

plt.xlabel('Haplotypes', fontsize=11)
plt.ylabel('Edit Distance', fontsize=11)
# axes[i].set_yticklabels([f'$10^{int(tick)}$' for tick in y_ticks])
plt.yticks(np.arange(0, 6, 1), [f'$10^{i}$' for i in range(6)])
plt.xticks(x, haplotypes, fontsize=11)
plt.grid(axis='y', linestyle='--', alpha=0.6, zorder=0)
plt.legend(fontsize=10, loc='upper center', ncol=5, bbox_to_anchor=(0.5, 1.20))
plt.tight_layout()

# Save figure
plt.savefig('increasing_edit.pdf', format='pdf', dpi=1200, bbox_inches='tight')
plt.close()

# Plot 2: Separate bar plots for runtime and memory usage
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8, 3))

# Plotting runtime
for i, cov in enumerate(coverage_levels):
    ax1.bar(x - (2-i)*width, runtimes[cov], width, label=f'{cov}', color=plt.get_cmap('tab10')(i), zorder=3)

ax1.set_xlabel('Haplotypes', fontsize=11)
ax1.set_ylabel('Runtime (hours)', fontsize=11)
ax1.set_xticks(x)
ax1.set_xticklabels(haplotypes, fontsize=11)
ax1.grid(axis='y', linestyle='--', alpha=0.6, zorder=0)

# Plotting memory usage
for i, cov in enumerate(coverage_levels):
    ax2.bar(x - (2-i)*width, memory_usages[cov], width, label=f'{cov}', color=plt.get_cmap('tab10')(i), zorder=3)

ax2.set_xlabel('Haplotypes', fontsize=11)
ax2.set_ylabel('Memory Usage (GB)', fontsize=11)
ax2.set_xticks(x)
ax2.set_xticklabels(haplotypes, fontsize=11)
ax2.grid(axis='y', linestyle='--', alpha=0.6, zorder=0)

# Adding a combined legend at the top
handles1, labels1 = ax1.get_legend_handles_labels()
fig.legend(handles1, labels1, loc='upper center', ncol=5, bbox_to_anchor=(0.5, 0.99), fontsize=10)

plt.tight_layout()
plt.subplots_adjust(top=0.85)

# Save figure
plt.savefig('increasing_performance.pdf', format='pdf', dpi=1200, bbox_inches='tight')
plt.close()