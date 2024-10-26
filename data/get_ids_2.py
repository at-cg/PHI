#!/usr/bin/env python3

import random
import sys

# Check if the arguments are provided
if len(sys.argv) != 3:
    print("Usage: python get_ids_2.py <ignore_genomes_string> <count>")
    sys.exit(1)

# Get the ignore_genomes string and count from command-line arguments
ignore_genomes = sys.argv[1].strip()
# print(f"ignore_genomes: {ignore_genomes}")  # Debugging print

try:
    count = int(sys.argv[2])
except ValueError:
    print("Error: <count> must be an integer.")
    sys.exit(1)

# Correctly extract the IDs from the ignore_genomes string
# Assumes the IDs are prefixed with '-R ' and separated by spaces
id_prefixes = [entry for entry in ignore_genomes.split() if not entry.startswith('-R')]
# print(f"Extracted IDs: {id_prefixes}")  # Debugging print

# Validate the count
if count < 1 or count > len(id_prefixes):
    print(f"Invalid count. Please enter a number between 1 and {len(id_prefixes)}.")
    sys.exit(1)

# Randomly select a subset of the IDs
random.seed()
selected_ids = random.sample(id_prefixes, count)

# Print the subsampled IDs in the same format as ignore_genomes
formatted_output = " ".join([f"-R {id_}" for id_ in selected_ids])
print(formatted_output)
