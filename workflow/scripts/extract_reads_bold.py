#!/usr/bin/env python3
"""
extract reads from bold data based on primers and format taxonomy.

this script replicates the logic from kim navarro-velez's r script for processing bold data. it reads a tsv file of bold data, finds primer-flanked amplicon regions, formats the taxonomy, and writes the output to fasta and tsv files for the next steps in the snakemake workflow.

usage:
    python extract_reads_bold.py \\
        --input-tsv input.tsv \\
        --output-fasta output.fasta \\
        --output-taxonomy output.taxonomy.tsv \\
        --forward-primer GGTCAACAAATCATAAAGATATTGG \\
        --reverse-primer GGWACTAATCAATTTCCAAATCC \\
        --min-length 100
"""

import sys
import pandas as pd
import re
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import argparse

def parse_arguments():
    """parse command-line arguments."""
    parser = argparse.ArgumentParser(description='extract amplicon reads from bold data and format taxonomy.')
    parser.add_argument('--input-tsv', required=True, help='input bold tsv file.')
    parser.add_argument('--output-fasta', required=True, help='output fasta file for extracted sequences.')
    parser.add_argument('--output-taxonomy', required=True, help='output tsv file for formatted taxonomy.')
    parser.add_argument('--forward-primer', required=True, help='forward primer sequence.')
    parser.add_argument('--reverse-primer', required=True, help='reverse primer sequence (5-3 orientation).')
    parser.add_argument('--locus', default='coi-5p', help='locus name for fasta headers (default: coi-5p).')
    parser.add_argument('--min-length', type=int, default=100, help='minimum length of sequence to keep after extraction.')
    return parser.parse_args()

def format_taxonomy(row):
    """
    format the taxonomy string to match the qiime-compatible format (e.g., k__;p__;c__;...).
    this replicates the `makefasta_function` from kim's original r script.
    """
    # we determine the kingdom based on the phylum. this is a simplification but
    # covers the major groups in kim's script.
    phylum = str(row.get('phylum_name', '')).lower()
    if phylum in ['ascomycota', 'basidiomycota', 'chytridiomycota', 'glomeromycota', 'myxomycota', 'zygomycota']:
        kingdom = "fungi"
    elif phylum in ['chlorarachniophyta', 'ciliophora', 'heterokontophyta', 'pyrrophycophyta']:
        kingdom = "protozoa"
    else:
        kingdom = "animalia"

    tax_levels = ['phylum', 'class', 'order', 'family', 'genus', 'species']
    tax_prefixes = {'kingdom': 'k__', 'phylum': 'p__', 'class': 'c__', 'order': 'o__', 'family': 'f__', 'genus': 'g__', 'species': 's__'}
    
    tax_string_parts = [f"{tax_prefixes['kingdom']}{kingdom}"]
    
    for level in tax_levels:
        level_name = f"{level}_name"
        value = str(row.get(level_name, ''))
        # ensure we have a valid, non-empty string before adding the prefix.
        if value and pd.notna(value):
            tax_string_parts.append(f"{tax_prefixes[level]}{value}")
        else:
            # if a rank is missing, we add an empty placeholder to maintain the hierarchy.
            tax_string_parts.append(f"{tax_prefixes[level]}")
            
    return ";".join(tax_string_parts)

def extract_amplicon_region(sequence, forward_primer, reverse_primer_original):
    """
    finds primers and extracts the dna sequence between them.
    handles ambiguous bases in primers (like 'w' in the reverse primer).
    """
    # the reverse primer is given in 5'-3' orientation, but for searching, we need its
    # reverse complement. biopython handles this, including ambiguous bases.
    fwd_primer_seq = Seq(forward_primer)
    rev_primer_seq = Seq(reverse_primer_original)
    rev_primer_rc = rev_primer_seq.reverse_complement()

    # create regular expression patterns that account for ambiguous dna codes.
    fwd_regex = str(fwd_primer_seq).replace('N', '.').replace('R', '[AG]').replace('Y', '[CT]').replace('S', '[GC]').replace('W', '[AT]').replace('K', '[GT]').replace('M', '[AC]').replace('B', '[CGT]').replace('D', '[AGT]').replace('H', '[ACT]').replace('V', '[ACG]')
    rev_regex = str(rev_primer_rc).replace('N', '.').replace('R', '[AG]').replace('Y', '[CT]').replace('S', '[GC]').replace('W', '[AT]').replace('K', '[GT]').replace('M', '[AC]').replace('B', '[CGT]').replace('D', '[AGT]').replace('H', '[ACT]').replace('V', '[ACG]')

    # search for the primers in the sequence string. we use re.search which finds
    # the first occurrence.
    fwd_match = re.search(fwd_regex, sequence, re.IGNORECASE)
    if not fwd_match:
        return None
    
    # we start searching for the reverse primer *after* the forward primer ends.
    rev_match = re.search(rev_regex, sequence, re.IGNORECASE)
    if not rev_match or rev_match.start() <= fwd_match.end():
        return None
    
    # if both are found, we extract the region between the start of the forward
    # primer and the end of the reverse primer.
    start_pos = fwd_match.start()
    end_pos = rev_match.end()
    
    return sequence[start_pos:end_pos]

def main():
    args = parse_arguments()
    
    print(f"reading bold data from {args.input_tsv}")
    
    try:
        df = pd.read_csv(args.input_tsv, sep='\t', low_memory=False)
        # we drop rows that don't have a sequence id or nucleotide sequence.
        df.dropna(subset=['sequenceid', 'nucleotides'], inplace=True)
    except Exception as e:
        print(f"error reading or processing tsv file: {e}", file=sys.stderr)
        sys.exit(1)
    
    print(f"processing {len(df)} records with sequences.")
    
    # these lists will store our successful results.
    fasta_records = []
    taxonomy_records = []
    
    for _, row in df.iterrows():
        sequence = row['nucleotides']
        # it's safer to convert sequenceid to string to avoid potential float issues.
        sequence_id = str(row['sequenceid'])
        
        # for each sequence, we try to extract the amplicon region.
        amplicon = extract_amplicon_region(sequence, args.forward_primer, args.reverse_primer)
        
        # we only proceed if an amplicon was found and it meets our minimum length requirement.
        if amplicon and len(amplicon) >= args.min_length:
            
            # create a biopython sequence record. this is a standard object
            # that holds a sequence and its metadata.
            seq_record = SeqRecord(
                Seq(amplicon),
                id=sequence_id,
                description="" # keep description empty for cleaner fasta file.
            )
            fasta_records.append(seq_record)
            
            # format the taxonomy for this record and store it.
            taxonomy = format_taxonomy(row)
            taxonomy_records.append(f"{sequence_id}\t{taxonomy}")
    
    print(f"successfully extracted {len(fasta_records)} sequences meeting criteria.")
    
    # write the fasta records to the output file.
    try:
        with open(args.output_fasta, 'w') as f_out:
            for record in fasta_records:
                f_out.write(f">{record.id}\n{record.seq}\n")
        print(f"wrote {len(fasta_records)} sequences to {args.output_fasta}")
    except IOError as e:
        print(f"error writing to fasta file {args.output_fasta}: {e}", file=sys.stderr)
        sys.exit(1)

    # write the taxonomy records to the output file.
    try:
        with open(args.output_taxonomy, 'w') as t_out:
            t_out.write("\n".join(taxonomy_records))
        print(f"wrote {len(taxonomy_records)} taxonomy entries to {args.output_taxonomy}")
    except IOError as e:
        print(f"error writing to taxonomy file {args.output_taxonomy}: {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main() 