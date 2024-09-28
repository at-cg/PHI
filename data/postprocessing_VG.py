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
        # Update hap_id with 'VG'
        hap_id = f'rec_hap_{read}_{cov}x_VG'
        log_file = f'{output_dir}{hap_id}.log'
        if os.path.exists(log_file):
            with open(log_file, 'r') as file:
                log_data = file.read()

            # Extract runtime (wall clock time)
            real_time_match = re.search(r'Elapsed \(wall clock\) time \(h:mm:ss or m:ss\):\s+([\d:]+)', log_data)
            real_time = real_time_match.group(1) if real_time_match else None

            # Convert real time to seconds if found
            if real_time:
                time_parts = real_time.split(":")
                if len(time_parts) == 3:  # h:mm:ss format
                    real_time = int(time_parts[0]) * 3600 + int(time_parts[1]) * 60 + float(time_parts[2])
                elif len(time_parts) == 2:  # m:ss format
                    real_time = int(time_parts[0]) * 60 + float(time_parts[1])

            # Extract peak RSS (look for "Maximum resident set size")
            peak_rss_match = re.search(r'Maximum resident set size \(kbytes\):\s+(\d+)', log_data)
            peak_rss = float(peak_rss_match.group(1)) / (1024 * 1024) if peak_rss_match else None  # Convert kB to GB

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
        hap_id = f'rec_hap_{read}_{cov}x_VG'
        table.append([hap_id] + data[read][hap_id])
    print(f"Read: {read}")
    print(tabulate(table, headers=['Haplotype', 'Real time(s)', 'Peak RSS(GB)', 'Edit distance']))
    print('\n\n')
    # Save the table to a file as text
    with open(f'{output_dir}stats_{read}.txt', 'w') as file:
        file.write(f"Read: {read}\n")
        file.write(tabulate(table, headers=['Haplotype', 'Real time(s)', 'Peak RSS(GB)', 'Edit distance']))
        file.write('\n\n')
