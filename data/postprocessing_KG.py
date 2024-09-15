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

# Function to convert elapsed time (format: M:SS.ss or H:MM:SS.ss) to seconds
def elapsed_time_to_seconds(elapsed):
    parts = elapsed.split(':')
    if len(parts) == 2:  # Format M:SS.ss
        minutes = int(parts[0])
        seconds = float(parts[1])
        total_seconds = minutes * 60 + seconds
    elif len(parts) == 3:  # Format H:MM:SS.ss
        hours = int(parts[0])
        minutes = int(parts[1])
        seconds = float(parts[2])
        total_seconds = hours * 3600 + minutes * 60 + seconds
    else:
        total_seconds = None
    return total_seconds

# Read log files and extract statistics
for read in reads:
    data[read] = {}
    for cov in coverage:
        hap_id = f'rec_hap_{read}_{cov}x_KG'
        log_file = f'{output_dir}{hap_id}.log'
        if os.path.exists(log_file):
            with open(log_file, 'r') as file:
                log_data = file.read()

            # Extract user, system, and elapsed time
            time_matches = re.findall(r'(\d+\.\d+)user\s+(\d+\.\d+)system\s+(\d+:\d+\.\d+)elapsed', log_data)
            total_real_time = 0
            if time_matches:
                for time_match in time_matches:
                    elapsed_time_str = time_match[2]  # Get elapsed time
                    real_time = elapsed_time_to_seconds(elapsed_time_str)  # Convert to seconds
                    total_real_time += real_time

            # Extract max resident set size (RSS)
            rss_matches = re.findall(r'(\d+)maxresident', log_data)
            max_rss = max([int(rss) for rss in rss_matches]) if rss_matches else None
            peak_rss_gb = max_rss / (1024 * 1024) if max_rss else None  # Convert to GB

            # Extract edit distance (if applicable)
            edit_distance_match = re.search(r'#0:\s+(\d+)\s+', log_data)
            edit_distance = int(edit_distance_match.group(1)) if edit_distance_match else None

            # Store the extracted data
            data[read][hap_id] = [total_real_time, peak_rss_gb, edit_distance]
        else:
            # Handle missing logs gracefully
            data[read][hap_id] = [None, None, None]

# Generate and print tables for each read and coverage level
for read in reads:
    table = []
    for cov in coverage:
        hap_id = f'rec_hap_{read}_{cov}x_KG'
        table.append([hap_id] + data[read][hap_id])
    print(f"Read: {read}")
    print(tabulate(table, headers=['Haplotype', 'Real time(s)', 'Peak RSS(GB)', 'Edit distance']))
    print('\n\n')
    # Save the table to a file as text
    with open(f'{output_dir}stats_{read}.txt', 'w') as file:
        file.write(f"Read: {read}\n")
        file.write(tabulate(table, headers=['Haplotype', 'Real time(s)', 'Peak RSS(GB)', 'Edit distance']))
        file.write('\n\n')
