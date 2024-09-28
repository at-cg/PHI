#!/usr/bin/env python3

import os
import subprocess
import argparse

def run_command(command):
    print(f"Running command: {command}")
    subprocess.run(command, shell=True, check=True)

def main(temp_gbz, read_file, tmp_dir, threads, output_fasta):
    # Create tmp_dir if it does not exist
    os.makedirs(tmp_dir, exist_ok=True)

    # Step 1: Create the index for the GBZ file
    run_command(f"vg index -j {tmp_dir}/temp.dist {temp_gbz}")

    # Step 2: Generate the reference index
    run_command(f"vg gbwt -p --num-threads {threads} -r {tmp_dir}/temp.ri -Z {temp_gbz}")

    # Step 3: Generate haplotypes
    run_command(f"vg haplotypes -v 2 -t {threads} -H {tmp_dir}/temp.hapl {temp_gbz}")

    # Step 4: Haplotype sampling using KMC
    run_command(f"kmc -k29 -m128 -okff -t{threads} -hp {read_file} {tmp_dir}/sample {tmp_dir}")

    # Step 5: Generate a single haplotype using the KMC database
    run_command(f"vg haplotypes -v 2 -t {threads} --num-haplotypes 1 "
                f"-i {tmp_dir}/temp.hapl -k {tmp_dir}/sample.kff -g {tmp_dir}/sample.gbz {temp_gbz}")

    # Step 6: Extract the paths to a FASTA file
    run_command(f"vg paths -x {tmp_dir}/sample.gbz -F -S recombination > {output_fasta}")

    # Step 7: Use seqtk to reverse complement the final FASTA file (overwrite the specified output file)
    run_command(f"seqtk seq -r {output_fasta} > temp.fa && mv temp.fa {output_fasta}")

    print(f"Process completed! The reverse strand FASTA file is saved as '{output_fasta}'.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process GBZ input and output a reverse strand FASTA file.")
    parser.add_argument("-g", "--gbz", type=str, required=True, help="Input GBZ file (e.g., temp.gbz)")
    parser.add_argument("-r", "--reads", type=str, required=True, help="Input read file (e.g., APD_10x.fastq)")
    parser.add_argument("-t", "--threads", type=int, default=16, help="Number of threads to use (default: 16)")
    parser.add_argument("-d", "--tmp-dir", type=str, required=True, help="Temporary directory for intermediate files")
    parser.add_argument("-o", "--output", type=str, default="sample.fa", help="Output FASTA file (default: sample.fa)")

    args = parser.parse_args()

    main(args.gbz, args.reads, args.tmp_dir, args.threads, args.output)
