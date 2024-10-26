#!/usr/bin/env python3

import os
import re
import subprocess
from tabulate import tabulate
from concurrent.futures import ThreadPoolExecutor, as_completed

# List of reads and coverage levels
reads = ["APD", "DBB", "MANN", "QBL", "SSTO"]
coverage = [0.1, 0.5, 1, 2, 5, 10, 15]

# Dictionary to hold extracted data
data = {}

# Ensure the output directory exists
output_dir = 'data/Rec_haps/'
os.makedirs(output_dir, exist_ok=True)

# Function to compute edit distance and alignment identity using edlib
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

        return edit_distance, alignment_identity
    except Exception as e:
        print(f"Error computing edlib metrics for {query_fasta}: {e}")
        return None, None

# Function to process log file for a given read and coverage level
def process_log_file(read, cov):
    hap_id = f'rec_hap_{read}_{cov}x_PG'
    log_file = f'{output_dir}{hap_id}.log'
    ground_truth_fasta = f'Ground_truth/{read}.fasta'
    query_fasta = f'{output_dir}{hap_id}.fa'

    if os.path.exists(log_file):
        with open(log_file, 'r') as file:
            log_data = file.read()

        # Extract runtime (total wallclock time)
        real_time_match = re.search(r'total wallclock time:\s+(\d+\.\d+)\s+sec', log_data)
        real_time = float(real_time_match.group(1)) if real_time_match else None

        # Extract peak RSS (look for "Max RSS")
        peak_rss_match = re.search(r'Max RSS:\s+(\d+\.\d+)\s+GB', log_data)
        peak_rss = float(peak_rss_match.group(1)) if peak_rss_match else None

        # Extract homozygous and heterozygous variant counts
        variants_match = re.search(r'Total variants:\s+(\d+), Homozygous variants:\s+(\d+), Heterozygous variants:\s+(\d+)', log_data)
        homo_variants = int(variants_match.group(2)) if variants_match else None
        het_variants = int(variants_match.group(3)) if variants_match else None

        # Compute edit distance and alignment identity in parallel
        with ThreadPoolExecutor() as executor:
            future = executor.submit(compute_edlib_metrics, ground_truth_fasta, query_fasta)
            edit_distance, alignment_identity = future.result()

        # Store the extracted data
        return [real_time, peak_rss, edit_distance, alignment_identity, homo_variants, het_variants]
    else:
        # Handle missing logs gracefully
        return [None, None, None, None, None, None]

# Use ThreadPoolExecutor to process log files in parallel
with ThreadPoolExecutor() as executor:
    futures = {}
    for read in reads:
        data[read] = {}
        for cov in coverage:
            futures[(read, cov)] = executor.submit(process_log_file, read, cov)

    # Collect results as they are completed
    for (read, cov), future in futures.items():
        hap_id = f'rec_hap_{read}_{cov}x_PG'
        data[read][hap_id] = future.result()

# Generate and print tables for each read and coverage level
for read in reads:
    table = []
    for cov in coverage:
        hap_id = f'rec_hap_{read}_{cov}x_PG'
        table.append([hap_id] + data[read][hap_id])
    print(f"Read: {read}")
    print(tabulate(table, headers=['Haplotype', 'Real time(s)', 'Peak RSS(GB)', 'Edit distance', 'Alignment Identity (%)', 'Homozygous Variants', 'Heterozygous Variants']))
    print('\n\n')
    # Save the table to a file as text
    with open(f'{output_dir}stats_{read}.txt', 'w') as file:
        file.write(f"Read: {read}\n")
        file.write(tabulate(table, headers=['Haplotype', 'Real time(s)', 'Peak RSS(GB)', 'Edit distance', 'Alignment Identity (%)', 'Homozygous Variants', 'Heterozygous Variants']))
        file.write('\n\n')