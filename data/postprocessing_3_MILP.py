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

# Read log files and extract statistics
for read in reads:
    data[read] = {}
    for cov in coverage:
        hap_id = f'rec_hap_{read}_{cov}x_2_milp'
        log_file = f'{output_dir}{hap_id}.log'
        if os.path.exists(log_file):
            with open(log_file, 'r') as file:
                log_data = file.read()

            # Extract recombination count
            recomb_cnt_match = re.search(r'Recombination count:\s+(\d+)', log_data)
            recomb_cnt = int(recomb_cnt_match.group(1)) if recomb_cnt_match else None

            # Extract other statistics
            real_time_match = re.search(r'Real time:\s+(\d+\.\d+)\s+sec', log_data)
            real_time = float(real_time_match.group(1)) if real_time_match else None

            peak_rss_match = re.search(r'Peak RSS:\s+(\d+\.\d+)\s+GB', log_data)
            peak_rss = float(peak_rss_match.group(1)) if peak_rss_match else None

            edit_distance_match = re.search(r'#0:\s+(\d+)\s+', log_data)
            edit_distance = int(edit_distance_match.group(1)) if edit_distance_match else None

            # Extract minimizers and ILP percentage
            minimizers_match = re.search(r'Indexed reads with spectrum size:\s+(\d+)', log_data)
            minimizers = int(minimizers_match.group(1)) if minimizers_match else None

            ilp_percentage_match = re.search(r'(\d+\.\d+)% Minimizers are in ILP', log_data)
            ilp_percentage = float(ilp_percentage_match.group(1)) if ilp_percentage_match else None

            # Extract filtered minimizers percentage
            filtered_minimizers_match = re.search(r'Filtered/Retained Minimizers:\s+(\d+\.\d+)/(\d+\.\d+)%', log_data)
            filtered_minimizers = float(filtered_minimizers_match.group(1)) if filtered_minimizers_match else None
            retained_minimizers = float(filtered_minimizers_match.group(2)) if filtered_minimizers_match else None

            # Store the extracted data
            data[read][hap_id] = [recomb_cnt, real_time, peak_rss, edit_distance, minimizers, ilp_percentage, filtered_minimizers]
        else:
            # Handle missing logs gracefully
            data[read][hap_id] = [None, None, None, None, None, None, None]

# Generate and print tables for each read and coverage level
for read in reads:
    table = []
    for cov in coverage:
        hap_id = f'rec_hap_{read}_{cov}x_2_milp'
        table.append([hap_id] + data[read][hap_id])
    print(f"Read: {read}")
    print(tabulate(table, headers=['Haplotype', 'Recombination Count', 'Real time(s)', 'Peak RSS(GB)', 'Edit distance', 'Minimizers (Reads)', '% Minimizers in ILP', '% Filtered Minimizers']))
    print('\n\n')
    # Save the table to a file as text
    with open(f'{output_dir}stats_{read}_3.txt', 'w') as file:
        file.write(f"Read: {read}\n")
        file.write(tabulate(table, headers=['Haplotype', 'Recombination Count', 'Real time(s)', 'Peak RSS(GB)', 'Edit distance', 'Minimizers (Reads)', '% Minimizers in ILP', '% Filtered Minimizers']))
        file.write('\n\n')