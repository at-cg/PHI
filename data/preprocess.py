#!/usr/bin/env python3

import os
import sys
import argparse
import multiprocessing
import time
import subprocess
import random

nthreads = 4

# pass threads from argument if there is no argument print the usage
if len(sys.argv) > 1:
    parser = argparse.ArgumentParser(description='Preprocess data')
    parser.add_argument('-t', '--threads', type=int, help='Number of threads to use for preprocessing')
    args = parser.parse_args()
    nthreads = args.threads
else:
    print("Usage: preprocess.py -t <number of threads>")
    sys.exit(0)


time_start = time.time()

# check if data folder exists, if not create it
if not os.path.exists("data"):
    os.system("mkdir -p data")
else:
    os.system("rm -rf data")
    os.system("mkdir -p data")

# download data
os.system("wget https://zenodo.org/records/6617246/files/MHC-61.agc?download=1 -O data/MHC-61.agc")
if not os.path.exists("data/hprc_haps"):
    os.system("mkdir -p data/hprc_haps")
else:
    os.system("rm -rf data/hprc_haps")
    os.system("mkdir -p data/hprc_haps")

# Extract data
os.system("agc getcol -o data/hprc_haps data/MHC-61.agc")

os.system("""awk '{gsub(/^>CHM13#0$/,">0"); print}' data/hprc_haps/MHC-CHM13.0.fa > temp.fasta && mv temp.fasta data/hprc_haps/MHC-CHM13.0.fa""")

# install minigraph cactus
os.system("sh install_cactus.sh")

# if js folder exists, remove it
if os.path.exists("js"):
    os.system("rm -rf js")

os.system("source data/cactus-bin-v2.9.0/venv-cactus-v2.9.0/bin/activate && cactus-pangenome ./js MHC.seqfile --outDir ./data/MHC-49_MC_out --outName MHC-49-MC  --reference CHM13.0 --vcf --maxCores 48 --indexCores 32 --mapCores 8 --batchSystem single_machine")
# move GFA files to data folder
os.system("sh chop_graph.sh") # chop the nodes into smaller pieces

# download reads and merge
# create a folder for reads
if not os.path.exists("data/reads"):
    os.system("mkdir -p data/reads")

# List of reads and their corresponding URLs
reads = {
    "APD": "https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR17272303/SRR17272303",
    "DBB": "https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR17272302/SRR17272302",
    "MANN": "https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR17272301/SRR17272301",
    "QBL": "https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR17272300/SRR17272300",
    "SSTO": "https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR17272299/SRR17272299"
}

# Function to download and process a read
def download_and_process_read(read_name, url):
    os.system(f"wget -O {read_name} {url} && fastq-dump --split-files --outdir data/reads {read_name} && cat data/reads/{read_name}_* > data/reads/{read_name}.fastq && rm data/reads/{read_name}_*")

# Download and process reads
pool = multiprocessing.Pool(processes=len(reads))
pool.starmap(download_and_process_read, reads.items())
pool.close()
pool.join()

reads = ["APD", "DBB", "MANN", "QBL", "SSTO"]

# downsample reads by covergae
coverage = [0.1, 0.5, 1, 2, 5, 10, 15]

# make a folder for downsampled reads
if not os.path.exists("data/reads_downsampled"):
    os.system("mkdir -p data/reads_downsampled")

mean_read_length = {}

for read in reads:
    # use seqkit to get the mean read length
    output = subprocess.check_output(f"source ~/.bashrc && seqkit stats data/reads/{read}.fastq", shell=True).decode('utf-8').strip()
    # Split the output into lines and extract the mean read length
    lines = output.splitlines()
    stats_line = lines[1]  # the second line contains the statistics for the fastq file
    avg_len = float(stats_line.split()[6])  # the 7th column contains the average read length
    mean_read_length[read] = avg_len

for read in reads:
    for cov in coverage:
        count_reads = int(cov * 5000000 / mean_read_length[read])
        if (cov == 15):
            count_reads = 1000000000 # use all the available reads
        seed = random.randint(1, 10000) + random.randint(10001, 100000) + random.randint(100001, 1000000) + random.randint(1000001, 10000000) + random.randint(10000001, 100000000)
        os.system(f"source ~/.bashrc && seqkit sample -s {seed} -n {count_reads} data/reads/{read}.fastq > data/reads_downsampled/{read}_{cov}x.fastq")

# Preprocess done
time_end = time.time()

print("Preprocessing done in", time_end - time_start, "seconds")
