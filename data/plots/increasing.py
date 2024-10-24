import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# Load the CSV file into a DataFrame
df_updated = pd.read_csv('increasing.csv')

# Add columns for 3M, 7M, 13M, 25M, and 49M with extracted data
coverage_levels = ["3H", "7H", "13H", "25H", "49H"]
legend_labels = ['3', '7', '13', '25', '49']  # The new labels for the legend

for cov in coverage_levels:
    df_updated[cov] = df_updated[cov].apply(lambda x: tuple(map(float, x.strip('()').split(', '))))

# Extract the relevant columns for plotting
haplotypes = df_updated["Read"].tolist()

edit_distances = {cov: df_updated[cov].apply(lambda x: x[2]).tolist() for cov in coverage_levels}
runtimes = {cov: df_updated[cov].apply(lambda x: x[0] / 3600).tolist() for cov in coverage_levels}  # Convert to hours
memory_usages = {cov: df_updated[cov].apply(lambda x: x[1]).tolist() for cov in coverage_levels}

x = np.arange(len(haplotypes))
width = 0.15  # Width of bars for the bar plot

# Create subplots for combining all three plots in a single figure
fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(11, 2))

# Plot 1: Edit Distance (Log Scale)
for i, cov in enumerate(coverage_levels):
    ax1.bar(x - (2-i)*width, np.log10(edit_distances[cov]), width, label=legend_labels[i], color=plt.get_cmap('tab10')(i), zorder=3)

ax1.set_xlabel('Haplotype', fontsize=11)
ax1.set_ylabel('Edit distance', fontsize=11)
ax1.set_yticks(np.arange(0, 6, 1))
ax1.set_yticklabels([f'$10^{i}$' for i in range(6)])
ax1.set_xticks(x)
ax1.set_xticklabels(haplotypes, fontsize=11, rotation=45)
ax1.grid(axis='y', linestyle='--', alpha=0.6, zorder=0)

# Plot 2: Runtime (hours)
for i, cov in enumerate(coverage_levels):
    ax2.bar(x - (2-i)*width, runtimes[cov], width, label=legend_labels[i], color=plt.get_cmap('tab10')(i), zorder=3)

ax2.set_xlabel('Haplotype', fontsize=11)
ax2.set_ylabel('Runtime (hours)', fontsize=11)
ax2.set_xticks(x)
ax2.set_xticklabels(haplotypes, fontsize=11, rotation=45)
ax2.grid(axis='y', linestyle='--', alpha=0.6, zorder=0)

# Plot 3: Memory Usage (GB)
for i, cov in enumerate(coverage_levels):
    ax3.bar(x - (2-i)*width, memory_usages[cov], width, label=legend_labels[i], color=plt.get_cmap('tab10')(i), zorder=3)

ax3.set_xlabel('Haplotype', fontsize=11)
ax3.set_ylabel('Memory usage (GB)', fontsize=11)
ax3.set_xticks(x)
ax3.set_xticklabels(haplotypes, fontsize=11, rotation=45)
ax3.grid(axis='y', linestyle='--', alpha=0.6, zorder=0)

# Add a single legend for all plots at the top with the updated label
handles, labels = ax1.get_legend_handles_labels()
# fig.legend(handles, labels, title='Number of reference haplotypes', loc='upper center', ncol=5, bbox_to_anchor=(0.5, 1.05), fontsize=10)
plt.text(0.31, 0.95, 'Number of reference haplotypes', fontsize=11, transform=fig.transFigure, ha='right')
fig.legend(handles, labels, loc='upper center', ncol=5, bbox_to_anchor=(0.53, 1.08), fontsize=12)

# Adjust layout to fit the legend and plots
plt.tight_layout(rect=[0, 0, 1, 0.95])
# plt.subplots_adjust(top=0.85)
plt.subplots_adjust(top=0.85, wspace=0.35, hspace=0.7)  # Added hspace for space between top and bottom plots

# Save the combined figure
plt.savefig('combined_figure_with_haplotypes.pdf', format='pdf', dpi=1200, bbox_inches='tight')
plt.close()