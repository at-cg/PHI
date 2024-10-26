#!/usr/bin/env python3

import random
import sys

# List of IDs as pairs (prefixes only)
id_prefixes = [
    "HG002", "HG00438", "HG005", "HG00621", "HG00741", "HG01106",
    "HG01109", "HG01123", "HG01258", "HG01358", "HG01361", "HG01891",
    "HG01928", "HG01952", "HG01978", "HG02080", "HG02257", "HG02486",
    "HG02559", "HG02622", "HG02717", "HG02886", "HG03540", "NA18906"
]

# Check if the count argument is provided
if len(sys.argv) != 2:
    print("Usage: python script_name.py <count>")
    sys.exit(1)

# Get the count of pairs from command-line arguments
try:
    count = int(sys.argv[1])
except ValueError:
    print("Error: <count> must be an integer.")
    sys.exit(1)

# Validate the count
if count < 1 or count > len(id_prefixes):
    print(f"Invalid count. Please enter a number between 1 and {len(id_prefixes)}.")
    sys.exit(1)

# Randomly select prefixes based on user input
random.seed()
selected_prefixes = random.sample(id_prefixes, count)

# Join the prefixes into a comma-separated string
ignore_genomes = "".join([f" -R {prefix}" for prefix in selected_prefixes])

# Output ignore_genomes string
# print("Generated ignore list for GFA:")
print(ignore_genomes)