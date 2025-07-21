#!/usr/bin/env python3
"""
script to gather bold data by systematically downloading bold data for comprehensive taxonomic coverage.
this script is an adaptation of kim navarro-velez's r script and is intended to be
run as part of the snakemake workflow.

usage: python gather_bold_data.py --output-dir output_dir --locus coi
"""

import sys
import os
import pandas as pd
import argparse
import time
from pathlib import Path
import bold

class BOLDDataGatherer:
    """comprehensive bold data gatherer following kim's approach."""
    
    def __init__(self, output_dir, locus="coi-5p"):
        self.output_dir = Path(output_dir)
        self.locus = locus.lower()
        
        # create output directories
        self.raw_dir = self.output_dir / "raw"
        self.processed_dir = self.output_dir / "processed"
        self.raw_dir.mkdir(parents=True, exist_ok=True)
        self.processed_dir.mkdir(parents=True, exist_ok=True)
        
    def query_bold(self, taxon, retries=3, delay=5):
        """
        query bold api for a specific taxon using the bold-py library.
        this is more robust than raw requests.
        """
        for attempt in range(retries):
            try:
                # the 'bold.api.get_data' function is used to retrieve specimen and sequence data
                # we specify the taxon name and indicate that we want data in tsv format.
                # 'include_data=True' ensures we get sequence data, not just specimen metadata.
                response = bold.api.get_data(taxon, format='tsv', include_data=True)
                if response:
                    return response
                else:
                    print(f"no data returned for {taxon}. it might be an empty taxon.")
                    return None
            except Exception as e:
                print(f"error querying bold for {taxon} (attempt {attempt + 1}/{retries}): {e}")
                if attempt < retries - 1:
                    print(f"retrying in {delay} seconds...")
                    time.sleep(delay)
                else:
                    print(f"failed to query bold for {taxon} after {retries} attempts.")
                    return None

    def save_bold_data(self, data, filename):
        """save bold data to a file in the 'raw' directory."""
        filepath = self.raw_dir / filename
        with open(filepath, 'w', encoding='utf-8') as f:
            f.write(data)
        print(f"saved {filename} with {len(data.splitlines())} lines")
        return filepath

    def process_and_combine_data(self, data_list, output_filename):
        """
        combines data from a list of tsv strings, filters by marker code,
        and saves to a single tsv file.
        """
        all_dfs = []
        for data in data_list:
            if data and data.strip():
                from io import StringIO
                try:
                    df = pd.read_csv(StringIO(data), sep='\t')
                    all_dfs.append(df)
                except pd.errors.ParserError as e:
                    print(f"warning: could not parse some data, skipping. error: {e}")
                except Exception as e:
                    print(f"an unexpected error occurred while reading data: {e}")

        if not all_dfs:
            print(f"no data to process for {output_filename}.")
            return

        combined_df = pd.concat(all_dfs, ignore_index=True)
        
        # kim's original script filters for the 'coi-5p' marker *after* downloading.
        # we replicate that logic here to ensure consistency.
        if 'markercode' in combined_df.columns:
            initial_count = len(combined_df)
            # case-insensitive filter for the specified locus (e.g., 'coi-5p')
            combined_df = combined_df[combined_df['markercode'].str.lower() == self.locus]
            filtered_count = len(combined_df)
            print(f"filtered by locus '{self.locus}': {initial_count} -> {filtered_count} records.")
        else:
            print("warning: 'markercode' column not found, cannot filter by locus.")

        # save the processed data to the 'processed' directory
        output_path = self.processed_dir / output_filename
        combined_df.to_csv(output_path, sep='\t', index=False)
        print(f"processed and saved data to {output_path}")

    def gather_taxa_list(self, taxa_list, output_filename):
        """gather data for a simple list of taxa."""
        print(f"\n=== gathering {output_filename.replace('.tsv', '').replace('_', ' ')} ===")
        all_data = []
        for taxon in taxa_list:
            print(f"querying {taxon}...")
            data = self.query_bold(taxon)
            if data:
                self.save_bold_data(data, f"{taxon.lower()}.tsv")
                all_data.append(data)
            time.sleep(2)  # be nice to bold api
        
        self.process_and_combine_data(all_data, output_filename)

    def gather_large_insect_orders_by_families(self):
        """
        gather data for the four largest insect orders, splitting by major families
        to avoid api limits, just like in kim's original script.
        """
        print("\n=== gathering large insect orders by families ===")
        
        # this dictionary defines the large orders and the major families within each
        # that we want to download individually.
        orders_and_families = {
            "coleoptera": ["carabidae", "chrysomelidae", "curculionidae", "staphylinidae"],
            "diptera": ["sciaridae", "cecidomyiidae", "chironomidae", "phoridae",
                        "muscidae", "culicidae", "ceratopogonidae", "tachinidae"],
            "hymenoptera": ["braconidae", "formicidae", "ichneumonidae", "platygastridae"],
            "lepidoptera": ["noctuidae", "erebidae", "sphingidae", "geometridae"]
        }
        
        for order, major_families in orders_and_families.items():
            self._gather_order_by_families(order, major_families)
    
    def _get_all_families_in_order(self, order_name):
        """
        helper function to get all families within a given insect order using bold's taxonomy api.
        this replicates the 'downstream' logic from kim's r script.
        """
        print(f"fetching all families for order: {order_name}...")
        try:
            # gets the taxonomic information for the specified order
            tax_info = bold.api.get_taxonomy(order_name, fuzzy=False)
            if tax_info and 'families' in tax_info:
                # we extract just the names of the families
                family_names = [f['taxon'] for f in tax_info['families']]
                print(f"found {len(family_names)} families in {order_name}.")
                return family_names
            else:
                print(f"could not retrieve families for {order_name}.")
                return []
        except Exception as e:
            print(f"error fetching taxonomy for {order_name}: {e}")
            return []

    def _gather_order_by_families(self, order_name, major_families):
        """
        gathers data for a single large order.
        1. fetches all families in the order.
        2. downloads specified 'major_families' individually.
        3. downloads all 'other' families.
        4. combines and processes them.
        """
        print(f"\n--- gathering {order_name.capitalize()} by families ---")
        
        all_families = self._get_all_families_in_order(order_name)
        if not all_families:
            print(f"skipping {order_name} as no families could be retrieved.")
            return

        # determine which families are 'other'
        major_families_set = set(f.lower() for f in major_families)
        other_families = [fam for fam in all_families if fam.lower() not in major_families_set]
        
        print(f"will download {len(major_families)} major families and {len(other_families)} other families.")

        all_data = []

        # gather major families
        for family in major_families:
            print(f"querying major family: {family}...")
            data = self.query_bold(family)
            if data:
                all_data.append(data)
                self.save_bold_data(data, f"{family.lower()}.tsv")
            time.sleep(2)

        # gather other families
        for family in other_families:
            print(f"querying other family: {family}...")
            data = self.query_bold(family)
            if data:
                all_data.append(data)
                # we don't save every single 'other' family to a file to avoid clutter
            time.sleep(2)
        
        self.process_and_combine_data(all_data, f"insects_{order_name.lower()}.tsv")
    
    def combine_all_processed_data(self):
        """
        combine all individual processed tsv files from the 'processed' directory
        into one final file.
        """
        print("\n=== combining all processed data files ===")
        
        all_files = list(self.processed_dir.glob("*.tsv"))
        
        if not all_files:
            print("no processed files to combine.")
            return None

        all_dfs = []
        for filepath in all_files:
            print(f"reading {filepath.name}...")
            try:
                df = pd.read_csv(filepath, sep='\t', low_memory=False)
                all_dfs.append(df)
            except Exception as e:
                print(f"error reading {filepath.name}, skipping. error: {e}")

        if not all_dfs:
            print("could not read any of the processed files.")
            return None

        final_df = pd.concat(all_dfs, ignore_index=True)
        
        # this is the final, comprehensive dataset ready for the next step in the pipeline.
        final_path = self.output_dir / "all_bold_data.tsv"
        final_df.to_csv(final_path, sep='\t', index=False)
        
        print(f"\ncomplete! all bold data gathered and combined into: {final_path}")
        print(f"total records in final combined file: {len(final_df)}")
        return final_path

    def run_complete_gathering(self):
        """run the complete bold data gathering process."""
        print("starting comprehensive bold data gathering...")
        
        # define the lists of taxa to gather, mirroring kim's script
        non_arthropod_chordate_animals = [
            "acanthocephala", "acoelomorpha", "annelida", "brachiopoda", 
            "bryozoa", "chaetognatha", "cnidaria", "ctenophora", "cycliophora",
            "echinodermata", "entoprocta", "gastrotricha", "gnathostomulida",
            "hemichordata", "kinorhyncha", "mollusca", "nematoda", 
            "nematomorpha", "nemertea", "onychophora", "phoronida", 
            "placozoa", "platyhelminthes", "porifera", "priapulida", 
            "rhombozoa", "rotifera", "sipuncula", "tardigrada", "xenacoelomorpha"
        ]
        self.gather_taxa_list(non_arthropod_chordate_animals, "non_arthropod_chordate_animals.tsv")

        fungi = ['ascomycota', 'basidiomycota', 'chytridiomycota', 'glomeromycota', 'myxomycota', 'zygomycota']
        self.gather_taxa_list(fungi, "fungi.tsv")

        protists = ['chlorarachniophyta', 'ciliophora', 'heterokontophyta', 'pyrrophycophyta']
        self.gather_taxa_list(protists, "protists.tsv")
        
        chordates = [
            "mammalia", "aves", "reptilia", "amphibia", "actinopterygii",
            "chondrichthyes", "cephalaspidomorphi", "myxini"
        ]
        self.gather_taxa_list(chordates, "chordates.tsv")

        arthropods_non_insects = [
            "arachnida", "malacostraca", "branchiopoda", "maxillopoda",
            "ostracoda", "diplopoda", "chilopoda", "symphyla", "pauropoda"
        ]
        self.gather_taxa_list(arthropods_non_insects, "arthropods_non_insects.tsv")
        
        insects_other_orders = [
            "archaeognatha", "zygentoma", "ephemeroptera", "odonata",
            "plecoptera", "dermaptera", "embioptera", "phasmatodea",
            "orthoptera", "mantodea", "blattodea", "isoptera",
            "psocoptera", "phthiraptera", "thysanoptera", "hemiptera",
            "neuroptera", "megaloptera", "raphidioptera", "strepsiptera",
            "mecoptera", "siphonaptera", "trichoptera"
        ]
        self.gather_taxa_list(insects_other_orders, "insects_other_orders.tsv")
        
        # handle the large insect orders by splitting them into families
        self.gather_large_insect_orders_by_families()
        
        # combine all the processed data into one final file
        self.combine_all_processed_data()

        print("\nscript finished.")
        return

def main():
    parser = argparse.ArgumentParser(description='comprehensive bold data gathering script.')
    parser.add_argument('--output-dir', required=True, help='directory to save the output files.')
    parser.add_argument('--locus', default='coi-5p', help='the marker locus to filter for (e.g., "coi-5p"). case-insensitive.')
    
    args = parser.parse_args()
    
    # a check to ensure the output directory from the snakemake rule is used
    # this script is designed to be called from snakemake, which passes the output directory
    # from the rule's 'params' block.
    if not os.path.exists(args.output_dir):
        print(f"creating output directory: {args.output_dir}")
        os.makedirs(args.output_dir)

    gatherer = BOLDDataGatherer(args.output_dir, args.locus)
    gatherer.run_complete_gathering()
    
    print(f"\nnext steps in the pipeline will use the combined data file located in:")
    print(f"{Path(args.output_dir) / 'all_bold_data.tsv'}")

if __name__ == "__main__":
    main() 