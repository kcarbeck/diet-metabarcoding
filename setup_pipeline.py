#!/usr/bin/env python3
"""
Setup script for the diet metabarcoding pipeline.
This script helps users configure their pipeline properly.
"""

import os
import sys
import yaml
from pathlib import Path

def check_config():
    """Check if the configuration file exists and is valid."""
    config_file = "config/orchards.yaml"
    if not os.path.exists(config_file):
        print(f"ERROR: Configuration file {config_file} not found!")
        return False
    
    try:
        with open(config_file, 'r') as f:
            config = yaml.safe_load(f)
        
        # Check required fields
        required_fields = ['project_name', 'samplesheet', 'metadata_tsv', 'locus']
        missing_fields = []
        for field in required_fields:
            if field not in config:
                missing_fields.append(field)
        
        if missing_fields:
            print(f"ERROR: Missing required configuration fields: {missing_fields}")
            return False
        
        print("✓ Configuration file is valid")
        return True
        
    except Exception as e:
        print(f"ERROR: Could not parse configuration file: {e}")
        return False

def check_data_files():
    """Check if required data files exist."""
    print("\nChecking data files...")
    
    # Check manifest file
    manifest_file = "data/metadata/orchards_manifest.tsv"
    if not os.path.exists(manifest_file):
        print(f"⚠️  WARNING: Manifest file {manifest_file} not found")
        print("   Please create this file with your sample data paths")
    else:
        print("✓ Manifest file exists")
    
    # Check metadata file
    metadata_file = "data/metadata/orchards_meta.tsv"
    if not os.path.exists(metadata_file):
        print(f"⚠️  WARNING: Metadata file {metadata_file} not found")
        print("   Please create this file with your sample metadata")
    else:
        print("✓ Metadata file exists")
    
    # Check if raw data directory has files
    raw_dir = "data/raw/orchards"
    if not os.path.exists(raw_dir):
        print(f"⚠️  WARNING: Raw data directory {raw_dir} not found")
        print("   Please create this directory and add your sequencing data")
    else:
        files = list(Path(raw_dir).glob("*.fastq*"))
        if not files:
            print(f"⚠️  WARNING: No fastq files found in {raw_dir}")
            print("   Please add your sequencing data files")
        else:
            print(f"✓ Found {len(files)} fastq files in raw data directory")

def check_environments():
    """Check if conda environments exist."""
    print("\nChecking conda environments...")
    
    env_files = [
        "envs/bold-pipeline.yml",
        "envs/qiime2-metagenome-2024.10.yml"
    ]
    
    for env_file in env_files:
        if os.path.exists(env_file):
            print(f"✓ Environment file exists: {env_file}")
        else:
            print(f"⚠️  WARNING: Environment file missing: {env_file}")

def check_scripts():
    """Check if required scripts exist."""
    print("\nChecking workflow scripts...")
    
    script_files = [
        "workflow/scripts/gather_bold_data.py",
        "workflow/scripts/validate_metadata.py",
        "workflow/scripts/clean_bold_fasta_and_taxonomy.py",
        "workflow/scripts/extract_alignment_region.py",
        "workflow/scripts/extract_reads_bold.py",
        "workflow/scripts/filter_by_geo.py",
        "workflow/scripts/rel_abund_plot.R"
    ]
    
    for script_file in script_files:
        if os.path.exists(script_file):
            print(f"✓ Script exists: {script_file}")
        else:
            print(f"⚠️  WARNING: Script missing: {script_file}")

def check_snakefile():
    """Check if main Snakefile exists."""
    print("\nChecking Snakefile...")
    
    if os.path.exists("Snakefile"):
        print("✓ Main Snakefile exists")
    else:
        print("⚠️  WARNING: Main Snakefile missing")

def main():
    """Run all checks."""
    print("Diet Metabarcoding Pipeline Setup Check")
    print("=" * 40)
    
    # Run all checks
    config_ok = check_config()
    check_data_files()
    check_environments()
    check_scripts()
    check_snakefile()
    
    print("\n" + "=" * 40)
    print("SETUP SUMMARY:")
    
    if config_ok:
        print("✓ Configuration is valid")
        print("\nNEXT STEPS:")
        print("1. Add your sequencing data to data/raw/orchards/")
        print("2. Update data/metadata/orchards_manifest.tsv with your sample paths")
        print("3. Update data/metadata/orchards_meta.tsv with your sample metadata")
        print("4. Set your email address in config/orchards.yaml (optional)")
        print("5. Install snakemake: conda install -c conda-forge -c bioconda snakemake")
        print("6. Run the pipeline: snakemake --use-conda")
    else:
        print("✗ Configuration has issues - please fix them before proceeding")
        sys.exit(1)

if __name__ == "__main__":
    main() 