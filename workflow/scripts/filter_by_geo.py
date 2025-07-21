#!/usr/bin/env python3
"""
filter bold data by geographic location.

this script filters a bold tsv data file to include only records from
specified geographic locations (countries). crucially, it retains any records
where geographic information is missing, as this is common in bold.

this script is designed to be called from the snakemake workflow.

usage:
    python filter_by_geo.py \\
        --input-file path/to/input.tsv \\
        --output-file path/to/output.tsv \\
        --countries "USA,Canada"
"""

import sys
import pandas as pd
import argparse
from pathlib import Path

def parse_arguments():
    """parse command-line arguments."""
    parser = argparse.ArgumentParser(description='filter bold data by geographic location.')
    parser.add_argument(
        '--input-file', 
        required=True, 
        type=Path,
        help='input bold tsv file.'
    )
    parser.add_argument(
        '--output-file', 
        required=True, 
        type=Path,
        help='output path for the filtered tsv file.'
    )
    # NOTE: --countries is no longer required so that users can choose
    # to filter solely by province/state or bounding box.
    parser.add_argument(
        '--countries',
        required=False,
        default='',
        type=str,
        help='optional: comma-separated list of countries to keep (e.g., "USA,Canada").'
    )

    parser.add_argument(
        '--states',
        required=False,
        default='',
        type=str,
        help='optional: comma-separated list of states/provinces to keep (e.g., "California,Ontario").'
    )

    parser.add_argument(
        '--bbox',
        required=False,
        default='',
        type=str,
        help='''optional: latitude_min,latitude_max,longitude_min,longitude_max bounding box.
                 example: "-15,15,170,190". records whose lat/lon are missing are retained.'''
    )
    return parser.parse_args()

def filter_by_geo(input_file, output_file, countries_str='', states_str='', bbox_str=''):
    """Filter a BOLD TSV file by geographic criteria.

    A record is *kept* unless it **explicitly** violates one of the provided
    filters.  Missing data never cause removal.

    1. Countries: keep if country is NA or in allow-list.
    2. States/Provinces: same logic on the `province_state` column.
    3. Bounding box: keep if latitude/longitude are NA or fall inside the box.

    All three criteria are combined with logical AND (i.e., a row must satisfy
    every filter that was supplied). If a filter is not provided it is ignored.
    """
    print(f"reading bold data from: {input_file}")
    
    # check if the input file exists
    if not input_file.exists():
        print(f"error: input file not found at {input_file}", file=sys.stderr)
        sys.exit(1)
        
    try:
        df = pd.read_csv(input_file, sep='\t', low_memory=False)
    except Exception as e:
        print(f"error reading tsv file: {e}", file=sys.stderr)
        sys.exit(1)
    
    original_count = len(df)
    print(f"original record count: {original_count}")

    # Build a boolean mask that starts as all True and progressively removes rows
    keep_mask = pd.Series(True, index=df.index)

    # ------------------------- COUNTRY FILTER -----------------------------
    if countries_str:
        if 'country' not in df.columns:
            print("warning: --countries supplied but no 'country' column found → ignoring country filter.")
        else:
            countries_to_keep = [c.strip().lower() for c in countries_str.split(',') if c.strip()]
            print(f"keeping records from countries: {countries_to_keep} (others with explicit country will be dropped)")

            keep_mask &= (df['country'].isnull() | df['country'].str.lower().isin(countries_to_keep))

    # ----------------------- STATE/PROVINCE FILTER -------------------------
    if states_str:
        if 'province_state' not in df.columns:
            print("warning: --states supplied but no 'province_state' column found → ignoring state filter.")
        else:
            states_to_keep = [s.strip().lower() for s in states_str.split(',') if s.strip()]
            print(f"keeping records from states/provinces: {states_to_keep}")

            keep_mask &= (df['province_state'].isnull() | df['province_state'].str.lower().isin(states_to_keep))

    # --------------------------- BBOX FILTER ------------------------------
    if bbox_str:
        try:
            lat_min, lat_max, lon_min, lon_max = [float(x) for x in bbox_str.split(',')]
            print(f"bounding box: lat {lat_min}→{lat_max}, lon {lon_min}→{lon_max}")

            # Choose column names flexibly
            lat_col = 'lat' if 'lat' in df.columns else 'latitude' if 'latitude' in df.columns else None
            lon_col = 'lon' if 'lon' in df.columns else 'longitude' if 'longitude' in df.columns else None

            if not lat_col or not lon_col:
                print("warning: --bbox supplied but lat/lon columns not found → ignoring bbox filter.")
            else:
                # Ensure numeric, coerce errors to NaN so we keep those rows
                lat_vals = pd.to_numeric(df[lat_col], errors='coerce')
                lon_vals = pd.to_numeric(df[lon_col], errors='coerce')

                within_bbox = (
                    lat_vals.isnull() | lon_vals.isnull() |
                    ((lat_vals >= lat_min) & (lat_vals <= lat_max) &
                     (lon_vals >= lon_min) & (lon_vals <= lon_max))
                )
                keep_mask &= within_bbox
        except ValueError:
            print("error: --bbox must contain exactly four comma-separated numbers (lat_min,lat_max,lon_min,lon_max). ignoring bbox filter.")

    # ----------------------------------------------------------------------
    filtered_df = df[keep_mask]

    filtered_count = len(filtered_df)
    removed_count = original_count - filtered_count

    print(f"filtered record count: {filtered_count}")
    print(f"removed {removed_count} records that did not meet the geographic filters.")
    
    # save the resulting dataframe to the specified output file.
    try:
        output_file.parent.mkdir(parents=True, exist_ok=True)
        filtered_df.to_csv(output_file, sep='\t', index=False)
        print(f"saved filtered data to: {output_file}")
    except Exception as e:
        print(f"error writing to output file: {e}", file=sys.stderr)
        sys.exit(1)

def main():
    """main execution function."""
    args = parse_arguments()
    # call the unified filter
    filter_by_geo(
        args.input_file,
        args.output_file,
        countries_str=args.countries,
        states_str=args.states,
        bbox_str=args.bbox,
    )

if __name__ == "__main__":
    main() 