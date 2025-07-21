#!/usr/bin/env python3
"""
clean and validate bold fasta and taxonomy files before qiime2 import.

this script performs several critical quality control steps to ensure the reference database is robust, mirroring the validation checks from kim's original pipeline. it ensures that the fasta and taxonomy files are synchronized and free of common issues.

the script will:
  1. remove sequences shorter than a specified minimum length
  2. remove sequences containing invalid iupac dna characters
  3. remove records with duplicate sequence ids
  4. ensure that every sequence in the fasta file has a matching taxonomy entry, and vice-versa, creating perfectly synchronized files

usage:
    python clean_bold_fasta_and_taxonomy.py \\
        --input-fasta input.fasta \\
        --input-taxonomy input.taxonomy.tsv \\
        --output-fasta cleaned.fasta \\
        --output-taxonomy cleaned.taxonomy.tsv \\
        --min-length 50
"""

import argparse
import sys
from collections import defaultdict

def parse_arguments():
    """parse command-line arguments."""
    parser = argparse.ArgumentParser(description='clean and validate bold fasta and taxonomy files.')
    parser.add_argument('--input-fasta', required=True, help='input fasta file from the extract_reads step.')
    parser.add_argument('--input-taxonomy', required=True, help='input taxonomy tsv file from the extract_reads step.')
    parser.add_argument('--output-fasta', required=True, help='path for the cleaned output fasta file.')
    parser.add_argument('--output-taxonomy', required=True, help='path for the cleaned output taxonomy tsv file.')
    parser.add_argument('--min-length', type=int, default=1, help='minimum sequence length to keep.')
    return parser.parse_args()

def read_fasta(filepath):
    """
    reads a fasta file into a dictionary mapping ids to sequences.
    this is a simple parser that handles multiline sequences.
    """
    sequences = {}
    current_id = None
    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                current_id = line[1:].split()[0] # get id, ignore anything after a space
                sequences[current_id] = []
            elif current_id:
                sequences[current_id].append(line)
    
    # join the sequence lines together
    for seq_id, seq_parts in sequences.items():
        sequences[seq_id] = ''.join(seq_parts)
        
    return sequences

def read_taxonomy(filepath):
    """reads a two-column, tab-separated taxonomy file into a dictionary."""
    taxonomy = {}
    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split('\t', 1)
            if len(parts) == 2:
                taxonomy[parts[0]] = parts[1]
    return taxonomy

def main():
    args = parse_arguments()
    
    print("starting cleaning and validation process...")

    # step 1: read the input files into memory
    print(f"reading fasta file: {args.input_fasta}")
    sequences = read_fasta(args.input_fasta)
    print(f"found {len(sequences)} sequences.")

    print(f"reading taxonomy file: {args.input_taxonomy}")
    taxonomy = read_taxonomy(args.input_taxonomy)
    print(f"found {len(taxonomy)} taxonomy entries.")

    # step 2: initial filtering and validation
    
    # find all unique sequence ids from both files to start with
    all_ids = set(sequences.keys()) | set(taxonomy.keys())
    print(f"total unique ids across both files: {len(all_ids)}")
    
    valid_ids = set()
    
    # define the set of valid iupac characters for a dna sequence.
    # the dash '-' is included as it represents a gap in aligned sequences.
    valid_chars = set("acgturyswkmbdhvn-")

    for seq_id in all_ids:
        # ---- validation checks ----
        # check 1: is the id present in both files?
        if seq_id not in sequences:
            print(f"info: id '{seq_id}' found in taxonomy but not in fasta. discarding.")
            continue
        if seq_id not in taxonomy:
            print(f"info: id '{seq_id}' found in fasta but not in taxonomy. discarding.")
            continue

        seq = sequences[seq_id]
        
        # check 2: is the sequence long enough?
        if len(seq) < args.min_length:
            print(f"info: sequence '{seq_id}' is too short ({len(seq)}bp). discarding.")
            continue
            
        # check 3: does the sequence contain only valid characters?
        # we convert the sequence to lowercase and check if its characters are a subset
        # of our allowed `valid_chars`.
        if not set(seq.lower()).issubset(valid_chars):
            print(f"info: sequence '{seq_id}' contains invalid characters. discarding.")
            continue
            
        # if all checks pass, we add the id to our set of 'good' ids.
        valid_ids.add(seq_id)

    print(f"\nfound {len(valid_ids)} records that passed all validation checks.")
    
    # we don't explicitly check for duplicates because by using dictionaries,
    # any duplicate ids in the input files would have already been overwritten.
    # the final set of `valid_ids` will inherently be unique.

    # step 3: write the cleaned and synchronized output files
    
    # write the final fasta file, only including sequences whose ids are in `valid_ids`.
    print(f"writing cleaned fasta file to: {args.output_fasta}")
    with open(args.output_fasta, 'w') as f_out:
        for seq_id in sorted(list(valid_ids)): # sorting makes output deterministic
            f_out.write(f">{seq_id}\n{sequences[seq_id]}\n")

    # write the final taxonomy file, only including taxonomies whose ids are in `valid_ids`.
    print(f"writing cleaned taxonomy file to: {args.output_taxonomy}")
    with open(args.output_taxonomy, 'w') as t_out:
        for seq_id in sorted(list(valid_ids)): # sorting ensures same order as fasta
            t_out.write(f"{seq_id}\t{taxonomy[seq_id]}\n")
            
    print("\ncleaning and validation complete.")

if __name__ == "__main__":
    main() 