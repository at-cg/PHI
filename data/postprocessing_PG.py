#!/usr/bin/env python3

import os
import re
from tabulate import tabulate

# List of reads and coverage levels
reads = ["APD", "DBB", "MANN", "QBL", "SSTO"]
coverage = [0.1, 0.5, 1, 2, 5, 10, 15]

# Dictionary to hold extracted data
data = {}

# Ensure the output directory exists
output_dir = 'data/Rec_haps/'
os.makedirs(output_dir, exist_ok=True)

# Source the bashrc file
os.system("source ~/.bashrc")

# Read log files and extract statistics
for read in reads:
    data[read] = {}
    for cov in coverage:
        hap_id = f'rec_hap_{read}_{cov}x_PG'
        log_file = f'{output_dir}{hap_id}.log'
        if os.path.exists(log_file):
            with open(log_file, 'r') as file:
                log_data = file.read()

            # Extract runtime (total wallclock time)
            real_time_match = re.search(r'total wallclock time:\s+(\d+\.\d+)\s+sec', log_data)
            real_time = float(real_time_match.group(1)) if real_time_match else None

            # Extract peak RSS (look for "Max RSS")
            peak_rss_match = re.search(r'Max RSS:\s+(\d+\.\d+)\s+GB', log_data)
            peak_rss = float(peak_rss_match.group(1)) if peak_rss_match else None

            # Extract edit distance
            edit_distance_match = re.search(r'#0:\s+(\d+)\s+', log_data)
            edit_distance = int(edit_distance_match.group(1)) if edit_distance_match else None

            # Store the extracted data
            data[read][hap_id] = [real_time, peak_rss, edit_distance]
        else:
            # Handle missing logs gracefully
            data[read][hap_id] = [None, None, None]

# Generate and print tables for each read and coverage level
for read in reads:
    table = []
    for cov in coverage:
        hap_id = f'rec_hap_{read}_{cov}x_PG'
        table.append([hap_id] + data[read][hap_id])
    print(f"Read: {read}")
    print(tabulate(table, headers=['Haplotype', 'Real time(s)', 'Peak RSS(GB)', 'Edit distance']))
    print('\n\n')
    # Save the table to a file as text
    with open(f'{output_dir}stats_{read}.txt', 'w') as file:
        file.write(f"Read: {read}\n")
        file.write(tabulate(table, headers=['Haplotype', 'Real time(s)', 'Peak RSS(GB)', 'Edit distance']))
        file.write('\n\n')
