import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import ast  # To safely evaluate string tuples

# Function to extract runtime and memory from the tuple string
def parse_tuple_string(tuple_string):
    try:
        return ast.literal_eval(tuple_string)
    except (ValueError, SyntaxError):
        return (0, 0, 0)

# Load data from CSV files
phi_data = pd.read_csv('PHI_MILP.csv')
phi_ilp_no_relax_data = pd.read_csv('PHI_ilp.csv')

# Define coverage labels (extended to include 0.1x, 0.5x, 10x)
coverage = ['0.1x', '0.5x', '1x', '2x', '5x', '10x', '15x']
total_coverage = len(coverage)

# Extract legend labels from the first column of the CSV
legend_labels = phi_data.iloc[:, 0]

# Parse the CSV columns to extract runtime and memory
def parse_data(data, coverage):
    num_coverage_in_data = data.shape[1] - 1  # Exclude the first column with labels
    limited_coverage = coverage[:num_coverage_in_data]
    return data.iloc[:, 1:num_coverage_in_data+1].applymap(parse_tuple_string), limited_coverage

phi_data_parsed, phi_coverage = parse_data(phi_data, coverage)
phi_ilp_no_relax_data_parsed, _ = parse_data(phi_ilp_no_relax_data, coverage)

# Extract runtime and memory separately
phi_runtime = phi_data_parsed.applymap(lambda x: x[0]) / 3600  # Convert runtime to hours
phi_ilp_no_relax_runtime = phi_ilp_no_relax_data_parsed.applymap(lambda x: x[0]) / 3600

phi_rss = phi_data_parsed.applymap(lambda x: x[1])
phi_ilp_no_relax_rss = phi_ilp_no_relax_data_parsed.applymap(lambda x: x[1])

# Calculate maximum y-limits based only on ilp data
max_runtime = max(phi_runtime.max().max(), phi_ilp_no_relax_runtime.max().max())
max_rss = max(phi_rss.max().max(), phi_ilp_no_relax_rss.max().max())

# Plot configuration
fig, axes = plt.subplots(1, 2, figsize=(12, 4))

# Data to plot and corresponding y-axis labels
plot_data = [(phi_runtime, 'Runtime (hours)'), (phi_ilp_no_relax_runtime, '')]

# Titles for each subplot
titles = ['(C) ILP (with relaxation)', '(D) ILP (no relaxation)']

# Width of the bars
bar_width = 0.15
x_indices = np.arange(len(phi_coverage))

# Loop through the axes and data
for ax, (data, ylabel), title in zip(axes.flatten(), plot_data, titles):
    for i, label in enumerate(legend_labels):
        ax.bar(x_indices + i * bar_width, data.iloc[i, :len(phi_coverage)], bar_width, label=label, zorder=3)
    
    ax.set_title(title, fontsize=15)
    ax.set_xlabel('Coverage', fontsize=15)
    ax.set_xticks(x_indices + bar_width * (len(legend_labels) / 2 - 0.5))
    ax.set_xticklabels(phi_coverage, fontsize=15, rotation=0)  # Rotate x-axis labels
    ax.grid(axis='y', linestyle='--', alpha=0.6, zorder=0)
    ax.tick_params(axis='y', labelsize=15)

    # Set y-axis limits based solely on ilp data
    ax.set_ylim(0, max_runtime * 1.1)  # Add 10% buffer to max
    yticks = [0, 1.5, 3, 4.5, 6, 7.5]
    ax.set_yticks(yticks)

    # Add y-axis label for both subplots
    ax.set_ylabel(ylabel, fontsize=15)

# Add the legend for the labels without the title
legend = fig.legend(legend_labels, loc='upper center', ncol=5, fontsize=15, bbox_to_anchor=(0.5, 1.09))
# Manually add a title on the left using plt.text, positioned relative to the figure
plt.text(0.14, 1.0, 'Haplotype', fontsize=15, transform=fig.transFigure, ha='center')

# Adjust layout
plt.tight_layout()
plt.subplots_adjust(top=0.85, wspace=0.2)

# compute average speed up for ilp vs Milp over all haplotypes
average_speedup = []
for j in range(len(phi_coverage)):
    for i in range(len(legend_labels)):
        average_speedup.append(phi_ilp_no_relax_runtime.iloc[i, j] / phi_runtime.iloc[i, j])
average_speedup = sum(average_speedup) / len(average_speedup)
print('Average speedup ILP vs MILP: %.2f' % average_speedup)

# Save the figure
plt.savefig('phi_ilp_vs_phi_milp.pdf', bbox_inches='tight', dpi=1200, format='pdf')
