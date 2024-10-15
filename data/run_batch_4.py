#!/usr/bin/env python3

import os
import sys
import argparse
import multiprocessing
import time

reads = ["APD", "DBB", "MANN"]
coverage = [0.1, 0.5, 1, 2, 5, 10, 15]

nbatches = 1  # just for testing
max_threads = multiprocessing.cpu_count()

time_start = time.time()

# pass threads from argument; if there is no argument, print the usage
if len(sys.argv) > 1:
    parser = argparse.ArgumentParser(description='Preprocess data')
    parser.add_argument('-b', '--batches', type=int, help='Number of batches to use for preprocessing')
    args = parser.parse_args()
    nbatches = args.batches
else:
    print("Usage: run.py -b <number of batches>")
    sys.exit(0)

running_threads = max_threads // nbatches
if running_threads < 1:  # corner case
    running_threads = 1

if (running_threads > 32):
    running_threads = 32

# if not exists, create a folder for Rec_haps
if not os.path.exists("data/Rec_haps"):
    os.system("mkdir -p data/Rec_haps")

def par_run_PHI(read_cov_pair):
    read, cov = read_cov_pair
    downsampled_read = f"data/reads_downsampled/{read}_{cov}x.fastq"
    output_file = f"data/Rec_haps/rec_hap_{read}_{cov}x_2.fa"
    log_file = f"data/Rec_haps/rec_hap_{read}_{cov}x_2.log"
    ground_truth = f"Ground_truth/{read}.fasta"
    
    cmd = (
        f"PHI -m0 -q1 -R100 -k31 -w25 -t{running_threads} "
        f"-g data/MHC_49-MC_30_2.gfa -r {downsampled_read} "
        f"-o {output_file} > {log_file} 2>&1 && "
        f"edlib-aligner {ground_truth} {output_file} >> {log_file}"
    )
    os.system(cmd)

# Create a list of (read, coverage) pairs
read_cov_pairs = [(read, cov) for read in reads for cov in coverage]

# use multiprocessing to map downsampled reads in parallel
pool = multiprocessing.Pool(nbatches)
pool.map(par_run_PHI, read_cov_pairs)
pool.close()
pool.join()

time_end = time.time()
print("Time taken: ", time_end - time_start)
