#!/usr/bin/env python3

import os
import sys
import subprocess
import argparse
import gzip
import shutil

def run_command(cmd, source_bashrc=False):
    if source_bashrc:
        cmd = f"source ~/.bashrc && {cmd}"
    try:
        subprocess.run(cmd, shell=True, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error while running command: {e}")
        sys.exit(1)

def decompress_file(input_file, output_file):
    if input_file.endswith(".gz"):
        with gzip.open(input_file, 'rb') as f_in:
            with open(output_file, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
    else:
        shutil.copyfile(input_file, output_file)

def main(vcf, ref):
    # Create temporary filenames for intermediate files
    temp_vcf = "temp_vcf.vcf.gz"
    temp_ref = "temp_ref.fasta"
    temp_vcf_bgzipped = "temp_vcf_bgzipped.vcf.gz"
    temp_xg = "temp.xg"
    temp_gbwt = "temp.gbwt"
    temp_paths_gbwt = "temp.paths.gbwt"
    temp_combined_gbwt = "temp.combined.gbwt"
    temp_gbz = "temp.gbz"
    temp_gfa = "temp.gfa"
    
    # Decompress input files if necessary
    decompress_file(vcf, temp_vcf)
    decompress_file(ref, temp_ref)

    # Update VCF and Reference, changing CHM13#0 -> REF#0
    run_command(f"awk 'BEGIN {{OFS=\"\\t\"}}; $1 == \"#CHROM\" {{found=1}} found && NR>1 && $1 !~ /^#/ {{sub($1, \"REF#0\", $1)}}1' {temp_vcf} | bgzip > {temp_vcf_bgzipped} && tabix -f -p vcf {temp_vcf_bgzipped}")
    run_command(f"awk '/^>/ {{$0=\">REF#0\"}} {{print}}' {temp_ref} > REF.fasta")
    run_command("samtools faidx REF.fasta")

    # Run the VG commands with sourcing bashrc
    run_command(f"vg construct -v {temp_vcf_bgzipped} -r REF.fasta -aS -p > {temp_xg}", source_bashrc=True)
    run_command(f"vg gbwt -x {temp_xg} -v {temp_vcf_bgzipped} -o {temp_gbwt}", source_bashrc=True)
    run_command(f"vg gbwt -x {temp_xg} -E -o {temp_paths_gbwt}", source_bashrc=True)
    run_command(f"vg gbwt -m {temp_gbwt} {temp_paths_gbwt} -o {temp_combined_gbwt}", source_bashrc=True)
    run_command(f"vg gbwt -x {temp_xg} {temp_combined_gbwt} --gbz-format -g {temp_gbz}", source_bashrc=True)
    run_command(f"gfa2gbwt -d temp -p -m 30", source_bashrc=True)

    # Output the final GFA to stdout
    with open(temp_gfa, 'r') as gfa_file:
        sys.stdout.write(gfa_file.read())

    # Clean up temporary files
    os.system(f"rm temp*")
    os.system("rm REF.fasta*")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate GFA from VCF and FASTA/FA files.")
    parser.add_argument("-v", "--vcf", required=True, help="Input VCF file (can be gzipped).")
    parser.add_argument("-r", "--ref", required=True, help="Input reference FASTA/FA file (can be gzipped).")

    args = parser.parse_args()

    main(args.vcf, args.ref)
