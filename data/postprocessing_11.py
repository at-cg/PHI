#!/usr/bin/env python3

import os
import re
from tabulate import tabulate

# List of reads and coverage levels
reads = ["APD", "DBB", "MANN", "QBL", "SSTO"]
coverage = [15]

# Dictionary to hold extracted data
data = {}

# Ensure the output directory exists
output_dir = 'data/Rec_haps/'
os.makedirs(output_dir, exist_ok=True)

# Function to compute edit distance and identity using edlib
def compute_edlib_metrics(ground_truth_fasta, query_fasta):
    try:
        # Run the commands: source bashrc, activate conda environment, and execute edlib_edits.py
        command = f'source ~/.bashrc && conda activate edlib && python3 edlib_edits.py {ground_truth_fasta} {query_fasta}'
        result = subprocess.run(command, shell=True, capture_output=True, text=True, executable='/bin/bash')

        # Extract the relevant values from the output
        edit_distance_match = re.search(r'Edit distance:\s+(\d+)', result.stdout)
        alignment_identity_match = re.search(r'Alignment identity:\s+(\d+\.\d+)%', result.stdout)

        edit_distance = int(edit_distance_match.group(1)) if edit_distance_match else None
        alignment_identity = float(alignment_identity_match.group(1)) if alignment_identity_match else None

        if edit_distance is None or alignment_identity is None:
            print(f"Error: Edit distance or Alignment identity not computed for {query_fasta}")
        else:
            print(f"Computed for {query_fasta}: Edit Distance = {edit_distance}, Alignment Identity = {alignment_identity}%")

        return edit_distance, alignment_identity
    except Exception as e:
        print(f"Error computing edlib metrics for {query_fasta}: {e}")
        return None, None

# Function to process each read and coverage level
def process_read_coverage(read, cov):
    hap_id = f'rec_hap_{read}_{cov}x_49_2'
    log_file = f'{output_dir}{hap_id}.log'
    ground_truth_fasta = f'Ground_truth/{read}.fasta'
    query_fasta = f'{output_dir}{hap_id}.fa'

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

        # Compute edit distance and alignment identity using edlib
        edit_distance, alignment_identity = compute_edlib_metrics(ground_truth_fasta, query_fasta)

        # Extract minimizers and ILP percentage
        minimizers_match = re.search(r'Indexed reads with spectrum size:\s+(\d+)', log_data)
        minimizers = int(minimizers_match.group(1)) if minimizers_match else None

        ilp_percentage_match = re.search(r'(\d+\.\d+)% Minimizers are in ILP', log_data)
        ilp_percentage = float(ilp_percentage_match.group(1)) if ilp_percentage_match else None

        # Extract filtered minimizers percentage
        filtered_minimizers_match = re.search(r'Filtered/Retained Minimizers:\s+(\d+\.\d+)/(\d+\.\d+)%', log_data)
        filtered_minimizers = float(filtered_minimizers_match.group(1)) if filtered_minimizers_match else None
        retained_minimizers = float(filtered_minimizers_match.group(2)) if filtered_minimizers_match else None

        # Return the data for this read and coverage
        return read, cov, [recomb_cnt, real_time, peak_rss, edit_distance, alignment_identity, minimizers, ilp_percentage, filtered_minimizers]
    else:
        # Handle missing logs gracefully
        return read, cov, [None, None, None, None, None, None, None, None]

# Use ThreadPoolExecutor to parallelize compute_edlib_metrics
with ThreadPoolExecutor() as executor:
    futures = []
    for read in reads:
        data[read] = {}
        for cov in coverage:
            futures.append(executor.submit(process_read_coverage, read, cov))

    # Collect results as they are completed
    for future in as_completed(futures):
        read, cov, result = future.result()
        hap_id = f'rec_hap_{read}_{cov}x_49_2'
        data[read][hap_id] = result

# Generate and print tables for each read and coverage level
for read in reads:
    table = []
    for cov in coverage:
        hap_id = f'rec_hap_{read}_{cov}x_49_2'
        table.append([hap_id] + data[read][hap_id])
    print(f"Read: {read}")
    print(tabulate(table, headers=['Haplotype', 'Recombination Count', 'Real time(s)', 'Peak RSS(GB)', 'Edit distance', 'Alignment Identity (%)', 'Minimizers (Reads)', '% Minimizers in ILP', '% Filtered Minimizers']))
    print('\n\n')
    # Save the table to a file as text
    with open(f'{output_dir}stats_{read}_49_2.txt', 'w') as file:
        file.write(f"Read: {read}\n")
        file.write(tabulate(table, headers=['Haplotype', 'Recombination Count', 'Real time(s)', 'Peak RSS(GB)', 'Edit distance', 'Alignment Identity (%)', 'Minimizers (Reads)', '% Minimizers in ILP', '% Filtered Minimizers']))
        file.write('\n\n')