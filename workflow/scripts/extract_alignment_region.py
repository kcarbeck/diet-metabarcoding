#!/usr/bin/env python3
"""
extract a specific region from an aligned fasta file.

this script uses the biopython library to reliably parse a fasta file and
extract a sub-sequence from each record based on start and end coordinates.
it is a key part of the alignment workflow, used to trim the full alignment
down to just the amplicon region of interest.

the coordinates are 0-based, and the end coordinate is exclusive, which is
the standard for python slicing. for example, to extract from position 61
to 467 (1-based, inclusive), you would use start=60 and end=467.

usage:
    python extract_alignment_region.py \\
        -i input_aligned.fasta \\
        -o output_trimmed.fasta \\
        -s 60 \\
        -e 467
"""

import argparse
import sys
from pathlib import Path
from Bio import SeqIO

def parse_arguments():
    """parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description='extract a specific region from aligned fasta sequences.'
    )
    parser.add_argument(
        '-i', '--input', 
        required=True, 
        type=Path,
        help='input fasta file with aligned sequences.'
    )
    parser.add_argument(
        '-o', '--output', 
        required=True, 
        type=Path,
        help='output fasta file for the trimmed regions.'
    )
    parser.add_argument(
        '-s', '--start', 
        type=int, 
        required=True, 
        help='start position for trimming (0-based, inclusive).'
    )
    parser.add_argument(
        '-e', '--end', 
        type=int, 
        required=True, 
        help='end position for trimming (0-based, exclusive).'
    )
    return parser.parse_args()

def extract_region(input_file, output_file, start_pos, end_pos):
    """
    extracts a region from each sequence in a fasta file and writes to a new file.
    
    args:
        input_file (pathlib.path): path to the input fasta file.
        output_file (pathlib.path): path to the output fasta file.
        start_pos (int): start coordinate of the slice (0-based).
        end_pos (int): end coordinate of the slice (exclusive).
    """
    print(f"extracting region from position {start_pos} to {end_pos-1} (inclusive).")
    print(f"  - input: {input_file}")
    print(f"  - output: {output_file}")
    
    # basic validation of the coordinates.
    if start_pos < 0:
        print("error: start position must be a non-negative number.", file=sys.stderr)
        sys.exit(1)
    if end_pos <= start_pos:
        print("error: end position must be greater than the start position.", file=sys.stderr)
        sys.exit(1)
    
    trimmed_records = []
    records_processed = 0
    
    # using biopython's seqio.parse is the standard and most robust way to read fasta files.
    # it correctly handles various fasta formats, including multiline sequences.
    for record in SeqIO.parse(input_file, "fasta"):
        # `record.seq` is a biopython seq object. we can slice it just like a python string.
        # the slice `[start_pos:end_pos]` extracts the desired region.
        trimmed_seq = record.seq[start_pos:end_pos]
        
        # we create a new sequence record with the same id but with the trimmed sequence.
        # setting the description to empty cleans up the fasta header.
        new_record = record
        new_record.seq = trimmed_seq
        new_record.description = ""
        trimmed_records.append(new_record)
        records_processed += 1
        
    if records_processed == 0:
        print("warning: no sequences were found in the input file.")
    
    # seqio.write takes a list of sequence records and writes them to a file.
    # this is much safer than writing line-by-line manually.
    try:
        SeqIO.write(trimmed_records, output_file, "fasta")
        print(f"successfully processed {records_processed} records.")
        print("region extraction complete.")
    except IOError as e:
        print(f"error: could not write to output file {output_file}. reason: {e}", file=sys.stderr)
        sys.exit(1)

def main():
    """main execution function."""
    args = parse_arguments()
    
    if not args.input.exists():
        print(f"error: input file '{args.input}' does not exist.", file=sys.stderr)
        sys.exit(1)
    
    try:
        # ensure the output directory exists before trying to write to it.
        args.output.parent.mkdir(parents=True, exist_ok=True)
        extract_region(args.input, args.output, args.start, args.end)
    except Exception as e:
        print(f"an unexpected error occurred: {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main() 