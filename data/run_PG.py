#!/usr/bin/env python3

import os
import sys
import argparse
import multiprocessing
import time

reads = ["APD", "DBB", "MANN", "QBL", "SSTO"]
coverage = [0.1, 0.5, 1, 2, 5, 10, 15]

nbatches = 1  # just for testing
max_threads = multiprocessing.cpu_count()

time_start = time.time()

# Pass threads from argument; if there is no argument, print the usage
if len(sys.argv) > 1:
    parser = argparse.ArgumentParser(description='Run PanGenie on reads')
    parser.add_argument('-b', '--batches', type=int, help='Number of batches to use for processing')
    args = parser.parse_args()
    nbatches = args.batches
else:
    print("Usage: run.py -b <number of batches>")
    sys.exit(0)

running_threads = max_threads // nbatches
if running_threads < 1:  # corner case
    running_threads = 1

if running_threads > 32:
    running_threads = 32

# If not exists, create a folder for Rec_haps
if not os.path.exists("data/Rec_haps"):
    os.system("mkdir -p data/Rec_haps")

def par_run_PanGenie(read_cov_pair):
    read, cov = read_cov_pair
    downsampled_read = f"data/reads_downsampled/{read}_{cov}x.fastq"
    output_prefix = f"data/Rec_haps/rec_hap_{read}_{cov}x_PG"
    ref_file = "data/hprc_haps/MHC-CHM13.0.fa"
    vcf_file = "data/MHC_49-MC.vcf"
    
    temp_dir = f"temp_{read}_{cov}_PG"
    os.makedirs(temp_dir, exist_ok=True)

    log_file = f"{output_prefix}.log"
    ground_truth = f"Ground_truth/{read}.fasta"

    # Source the bashrc file, run PanGenie, and append edit distance to the log
    cmd = (
        f"source ~/.bashrc && "
        f"PanGenie -i {downsampled_read} -r {ref_file} -v {vcf_file} -o {temp_dir}/{read}_PG -t {running_threads} > {log_file} 2>&1 && "
        f"bgzip {temp_dir}/{read}_PG_genotyping.vcf && "
        f"tabix -p vcf {temp_dir}/{read}_PG_genotyping.vcf.gz && "
        f"total=$(bcftools view -H {temp_dir}/{read}_PG_genotyping.vcf.gz | wc -l); homo=$(bcftools view -i 'GT=\"hom\"' {temp_dir}/{read}_PG_genotyping.vcf.gz | wc -l); het=$((total - homo)); echo \"Total variants: $total, Homozygous variants: $homo, Heterozygous variants: $het\" >> {log_file} && "
        f"bcftools view -i 'GT=\"hom\"' {temp_dir}/{read}_PG_genotyping.vcf.gz | bgzip > {temp_dir}/{read}_PG_genotyping_no_het.vcf.gz && "
        f"tabix -p vcf {temp_dir}/{read}_PG_genotyping_no_het.vcf.gz && "
        f"bcftools consensus -f {ref_file} -o {output_prefix}.fa {temp_dir}/{read}_PG_genotyping_no_het.vcf.gz && "
        f"edlib-aligner {ground_truth} {output_prefix}.fa >> {log_file}"
    )


    os.system(cmd)

    # Clean up temporary directory
    os.system(f"rm -rf {temp_dir}")

# Create a list of (read, coverage) pairs
read_cov_pairs = [(read, cov) for read in reads for cov in coverage]

# Use multiprocessing to run PanGenie on different read-coverage pairs in parallel
pool = multiprocessing.Pool(nbatches)
pool.map(par_run_PanGenie, read_cov_pairs)
pool.close()
pool.join()

time_end = time.time()
print("Time taken: ", time_end - time_start)