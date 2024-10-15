#!/usr/bin/env python3

import edlib
import argparse
from Bio import SeqIO

def read_fasta_sequence(file):
    """Reads the first sequence from a FASTA file."""
    with open(file, "r") as fasta_file:
        # Use SeqIO to parse the FASTA file and return the sequence of the first record
        return str(next(SeqIO.parse(fasta_file, "fasta")).seq)

def compute_alignment_identity(seq1, seq2):
    # Perform alignment using edlib
    result = edlib.align(seq1, seq2, task="path")

    # Get the edit distance from the result
    edit_distance = result['editDistance']

    # Compute alignment length
    alignment_length = len(seq1)  # This can be len(seq1) or len(seq2), both sequences are aligned

    # Compute alignment identity
    if alignment_length > 0:
        alignment_identity = (alignment_length - edit_distance) * 100 / alignment_length
    else:
        alignment_identity = 0  # Avoid division by zero if alignment length is 0

    return edit_distance, alignment_identity

def main():
    # Set up argument parser
    parser = argparse.ArgumentParser(description="Compute edit distance and alignment identity between two FASTA sequences.")

    # Define positional arguments for query and reference FASTA files
    parser.add_argument("query_fasta", help="Path to the query FASTA file")
    parser.add_argument("reference_fasta", help="Path to the reference FASTA file")

    # Parse the arguments
    args = parser.parse_args()

    # Read sequences from the two FASTA files
    seq1 = read_fasta_sequence(args.query_fasta)
    seq2 = read_fasta_sequence(args.reference_fasta)

    # Get edit distance and alignment identity
    edit_distance, alignment_identity = compute_alignment_identity(seq1, seq2)

    # Print the results
    print(f"Edit distance: {edit_distance}")
    print(f"Alignment identity: {alignment_identity:.2f}%")

if __name__ == "__main__":
    main()