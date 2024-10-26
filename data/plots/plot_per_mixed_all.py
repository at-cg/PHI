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
phi_data = pd.read_csv('PHI_MIQP.csv')
phi_ilp_data = pd.read_csv('PHI_MILP.csv')
phi_ilp_no_relax_data = pd.read_csv('PHI_ilp.csv')
phi_iqp_no_relax_data = pd.read_csv('PHI.csv')

# Define coverage labels (extended to include 0.1x, 0.5x, 10x)
all_coverage = ['0.1x', '0.5x', '1x', '2x', '5x', '10x', '15x']
# coverage = ['1x', '2x', '5x', '10x', '15x']
coverage = ['0.1x', '0.5x', '1x', '2x', '5x', '10x', '15x']
total_coverage = len(all_coverage)
offset = (total_coverage - len(coverage)) + 1  # Offset to include 15x

# Extract legend labels from the first column of the CSV
legend_labels = phi_data.iloc[:, 0]

# Parse the CSV columns to extract runtime and memory
def parse_data(data, coverage):
    # Ensure there is data for each coverage level, otherwise limit the array
    num_coverage_in_data = data.shape[1] - 1  # Exclude the first column with labels
    limited_coverage = coverage[:num_coverage_in_data]
    return data.iloc[:, offset:num_coverage_in_data+1].applymap(parse_tuple_string), limited_coverage

phi_data_parsed, phi_coverage = parse_data(phi_data, coverage)
phi_ilp_data_parsed, _ = parse_data(phi_ilp_data, coverage)
phi_ilp_no_relax_data_parsed, _ = parse_data(phi_ilp_no_relax_data, coverage)
phi_iqp_no_relax_data_parsed, _ = parse_data(phi_iqp_no_relax_data, coverage)

# Extract runtime and memory separately
phi_runtime = phi_data_parsed.applymap(lambda x: x[0]) / 3600  # Convert runtime to hours
phi_ilp_runtime = phi_ilp_data_parsed.applymap(lambda x: x[0]) / 3600  # Convert runtime to hours
phi_ilp_no_relax_runtime = phi_ilp_no_relax_data_parsed.applymap(lambda x: x[0]) / 3600
phi_iqp_no_relax_runtime = phi_iqp_no_relax_data_parsed.applymap(lambda x: x[0]) / 3600

phi_rss = phi_data_parsed.applymap(lambda x: x[1])
phi_ilp_rss = phi_ilp_data_parsed.applymap(lambda x: x[1])
phi_ilp_no_relax_rss = phi_ilp_no_relax_data_parsed.applymap(lambda x: x[1])
phi_iqp_no_relax_rss = phi_iqp_no_relax_data_parsed.applymap(lambda x: x[1])

# Calculate maximum y-limits for side-by-side plots
max_runtime = max(phi_runtime.max().max(), phi_ilp_runtime.max().max(), 
                  phi_ilp_no_relax_runtime.max().max(), phi_iqp_no_relax_runtime.max().max())
max_rss = max(phi_rss.max().max(), phi_ilp_rss.max().max(), 
              phi_ilp_no_relax_rss.max().max(), phi_iqp_no_relax_rss.max().max())

# Plot configuration
fig, axes = plt.subplots(2, 4, figsize=(12, 5.5))

# Data to plot and corresponding y-axis labels
plot_data = [(phi_ilp_runtime, 'Runtime (hours)'), (phi_ilp_no_relax_runtime, None),
             (phi_runtime, None), (phi_iqp_no_relax_runtime, None),
             (phi_ilp_rss, 'Memory usage (GB)'), (phi_ilp_no_relax_rss, None), 
             (phi_rss, None), (phi_iqp_no_relax_rss, None)]

# Titles for each subplot
titles = ['(A) ILP (with relaxation)', '(B) ILP (no relaxation)', 
          '(C) IQP (with relaxation)', '(D) IQP (no relaxation)', 
          '(E) ILP (with relaxation)', '(F) ILP (no relaxation)', 
          '(G) IQP (with relaxation)', '(H) IQP (no relaxation)']

# Width of the bars
bar_width = 0.15
x_indices = np.arange(len(phi_coverage))

# Loop through the axes and data
for ax, (data, ylabel), title in zip(axes.flatten(), plot_data, titles):
    for i, label in enumerate(legend_labels):
        # Plot only for the available coverage levels
        ax.bar(x_indices + i * bar_width, data.iloc[i, :len(phi_coverage)], bar_width, label=label, zorder=3)
    
    ax.set_title(title, fontsize=13)
    ax.set_xlabel('Coverage', fontsize=13)
    ax.set_xticks(x_indices + bar_width * (len(legend_labels) / 2 - 0.5))
    ax.set_xticklabels(phi_coverage, fontsize=13, rotation=90)  # Rotate x-axis labels
    ax.grid(axis='y', linestyle='--', alpha=0.6, zorder=0)
    ax.tick_params(axis='y', labelsize=13)

    # Set y-axis limits based on max values
    if '(A)' in title or '(B)' in title or '(C)' in title or '(D)' in title:
        ax.set_ylim(0, max_runtime * 1.1)  # Add 10% buffer to max
        # ax.set_yticks(np.linspace(0, max_runtime, 5))  # Add 5 ticks
        yticks = [0, 1.5, 3, 4.5, 6, 7.5]  # round to 1 decimal places
        ax.set_yticks(yticks)
    else:
        ax.set_ylim(0, max_rss * 1.1)
        yticks = [0, 50, 75, 100, 125]
        ax.set_yticks(yticks)
    
    # Only add y-axis labels to left plots
    if ylabel is not None:
        ax.set_ylabel(ylabel, fontsize=13)

# Add the legend for the labels without the title
legend = fig.legend(legend_labels, loc='upper center', ncol=5, fontsize=13, bbox_to_anchor=(0.5, 1.0))

# Manually add a title on the left using plt.text, positioned relative to the figure
plt.text(0.19, 0.95, 'Haplotype', fontsize=13, transform=fig.transFigure, ha='center')

# Adjust layout
plt.tight_layout()
plt.subplots_adjust(top=0.85, wspace=0.3, hspace=0.9)  # Added hspace for space between top and bottom plots

# Save the figure
plt.savefig('phi_vs_phi_miqp_milp_ilp_no_relax_all.pdf', bbox_inches='tight', dpi=1200, format='pdf')