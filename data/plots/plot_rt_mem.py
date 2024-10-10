import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import ast

# Function to safely evaluate string tuples from CSV
def parse_tuple_string(tuple_string):
    try:
        return ast.literal_eval(tuple_string)
    except (ValueError, SyntaxError):
        return (0, 0, 0)

# Load data from CSV files
vg_data = pd.read_csv('VG.csv')
pangenie_data = pd.read_csv('PanGenie.csv')
kage_data = pd.read_csv('KAGE.csv')

# Parsing the CSV columns to extract runtime (x) and memory usage (y) from tuples
def parse_data(data):
    return data.applymap(parse_tuple_string)

vg_parsed = parse_data(vg_data.iloc[:, 1:])
pangenie_parsed = parse_data(pangenie_data.iloc[:, 1:])
kage_parsed = parse_data(kage_data.iloc[:, 1:])

reads = vg_data['Read']  # Read labels (APD, DBB, MANN, QBL, SSTO)

# Extracting runtime (in minutes) and memory usage (GB)
vg_runtime = vg_parsed.applymap(lambda x: x[0] / 60)  # Convert runtime to minutes (x)
pangenie_runtime = pangenie_parsed.applymap(lambda x: x[0] / 60)  # Convert runtime to minutes (x)
kage_runtime = kage_parsed.applymap(lambda x: x[0] / 60)  # Convert runtime to minutes (x)

vg_memory = vg_parsed.applymap(lambda x: x[1])  # Memory usage in GB (y)
pangenie_memory = pangenie_parsed.applymap(lambda x: x[1])  # Memory usage in GB (y)
kage_memory = kage_parsed.applymap(lambda x: x[1])  # Memory usage in GB (y)

# Customizing the style and font
plt.rc('font', family='sans-serif') 
plt.rc('font', serif='Helvetica Neue') 
plt.rc('text', usetex='false') 
plt.rcParams.update({'font.size': 22})

# Increase font size for Seaborn
sns.set(font_scale=1.2)

# Create a custom color palette with transparency by modifying the alpha
base_palette = sns.color_palette("Set1", len(reads))
custom_palette = [(r, g, b, 0.9) for r, g, b in base_palette]  # Set transparency (alpha = 0.9)

# Create 6 box plots (3 for runtime and 3 for memory)
fig, axes = plt.subplots(2, 3, figsize=(10, 6))

# Plot configurations
datasets = ['VG', 'PanGenie', 'KAGE']
runtime_data = [vg_runtime, pangenie_runtime, kage_runtime]
memory_data = [vg_memory, pangenie_memory, kage_memory]

# Titles for the plots with (A), (B), etc.
titles = ['(A) VG', '(B) PanGenie', '(C) KAGE', 
          '(D) VG', '(E) PanGenie', '(F) KAGE']

# Apply Seaborn whitegrid style
sns.set(style="whitegrid")

# Runtime plots (top row)
for i, (ax, data, title) in enumerate(zip(axes[0], runtime_data, titles[:3])):
    sns.boxplot(data=data.values.T, ax=ax, palette=custom_palette, zorder=3)  # Removed `alpha` parameter
    ax.set_title(f'{title} runtime', fontsize=12)
    ax.set_xticklabels(reads, rotation=45, fontsize=12)  # Rotate x-axis labels
    ax.set_xlabel('MHC haplotypes', fontsize=12)  # Add x-axis label
    if i == 0:
        ax.set_ylabel('Runtime (minutes)', fontsize=12)  # Only set y-label for the first plot
    ax.set_ylim(0, max(vg_runtime.max().max(), pangenie_runtime.max().max(), kage_runtime.max().max()) * 1.1)  # Same y-scale
    ax.grid(True, color='white', linestyle='-', linewidth=1, zorder=1)  # Grid behind box plots

# Memory usage plots (bottom row)
for i, (ax, data, title) in enumerate(zip(axes[1], memory_data, titles[3:])):
    sns.boxplot(data=data.values.T, ax=ax, palette=custom_palette, zorder=3)  # Removed `alpha` parameter
    ax.set_title(f'{title} memory usage', fontsize=12)
    ax.set_xticklabels(reads, rotation=45, fontsize=12)
    ax.set_xlabel('MHC haplotypes', fontsize=12)  # Add x-axis label
    if i == 0:
        ax.set_ylabel('Memory usage (GB)', fontsize=12)  # Only set y-label for the first plot
    ax.set_ylim(0, max(vg_memory.max().max(), pangenie_memory.max().max(), kage_memory.max().max()) * 1.1)  # Same y-scale
    ax.grid(True, color='white', linestyle='-', linewidth=1, zorder=1)  # Grid behind box plots

# Adjust layout
plt.tight_layout()

# Save the figure
plt.savefig('runtime_memory_boxplots_sns.pdf', bbox_inches='tight', dpi=1200, format='pdf')