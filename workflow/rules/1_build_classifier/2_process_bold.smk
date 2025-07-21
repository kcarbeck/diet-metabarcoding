# File: workflow/rules/1_build_classifier/2_process_bold.smk
# Description: Rules for filtering, extracting, and cleaning BOLD data.

# This file contains the rules that correspond to the initial data processing
# steps in `Kims_classifier_pipeline.rmd`. This includes optional geographic
# filtering, extraction of amplicon regions based on primers, and a critical
# cleaning and validation step to ensure data integrity before importing to QIIME2.

# ---------------------------------------------------------------------
# Rule 2.1: Filter BOLD data by geography (optional)
# ---------------------------------------------------------------------
rule bold_geo_filter:
    """
    Filter BOLD data by geographic location.
    
    This rule applies geographic filters based on the settings in the project's
    config file. It can filter by country, state/province, or a bounding box.
    If no geographic filters are specified in the config, this rule will simply
    copy the input data to the output path, allowing the pipeline to proceed seamlessly.
    """
    input:
        # Depends on the raw data from the gathering step.
        raw_data = f"{workdir}/bold/raw_bold_data.tsv"
    output:
        # The output is a new TSV file with the filtered data.
        filtered_data = f"{workdir}/bold/bold_filtered.tsv"
    params:
        # Geographic filter parameters are pulled from the config.
        countries = lambda w: ','.join(config.get("geo_filter", {}).get("countries", [])),
        states = lambda w: ','.join(config.get("geo_filter", {}).get("states", [])),
        bbox = lambda w: ','.join(map(str, config.get("geo_filter", {}).get("bbox", [])))
    log:
        f"{workdir}/logs/bold_geo_filter.log"
    conda:
        "../../envs/bold-pipeline.yml"
    shell:
        """
        set -euo pipefail
        
        # If no geographic filters are specified, just copy the input to the output.
        if [ -z "{params.countries}" ] && [ -z "{params.states}" ] && [ -z "{params.bbox}" ]; then
            echo "No geographic filter specified. Copying raw data to filtered path."
            cp {input.raw_data} {output.filtered_data}
        else
            echo "Filtering BOLD data by geography..."
            cmd="python workflow/scripts/filter_by_geo.py --input-file {input.raw_data} --output-file {output.filtered_data}"
            
            # Add filters to the command if they are specified in the config.
            if [ -n "{params.countries}" ]; then
                cmd="$cmd --countries '{params.countries}'"
            fi
            if [ -n "{params.states}" ]; then
                cmd="$cmd --states '{params.states}'"
            fi
            if [ -n "{params.bbox}" ]; then
                cmd="$cmd --bbox '{params.bbox}'"
            fi
            
            echo "Command: $cmd"
            $cmd 2>&1 | tee {log}
        fi
        """

# ---------------------------------------------------------------------
# Rule 2.2: Extract reads and format taxonomy
# ---------------------------------------------------------------------
rule bold_extract_reads:
    """
    Extract amplicon regions and format taxonomy.
    
    This rule uses primer sequences defined in the config to extract the target
    amplicon region from the BOLD sequences. It also formats the taxonomic
    information into a QIIME2-compatible string, replicating the logic from
    Kim's R script.
    """
    input:
        # Depends on the geographically filtered (or unfiltered) data.
        tsv = f"{workdir}/bold/bold_filtered.tsv"
    output:
        # Produces a FASTA file of sequences and a TSV of taxonomies.
        fasta = f"{workdir}/bold/bold_trimmed_seqs.fasta",
        taxonomy = f"{workdir}/bold/bold_trimmed_taxonomy.tsv"
    params:
        f_primer = forward_primer,
        r_primer = reverse_primer,
        locus = locus,
        min_length = bold_min_length
    log:
        f"{workdir}/logs/bold_extract_reads.log"
    conda:
        "../../envs/bold-pipeline.yml"
    shell:
        """
        set -euo pipefail
        
        echo "Extracting amplicon sequences from BOLD data..."
        python workflow/scripts/extract_reads_bold.py \\
            --input-tsv {input.tsv} \\
            --output-fasta {output.fasta} \\
            --output-taxonomy {output.taxonomy} \\
            --forward-primer "{params.f_primer}" \\
            --reverse-primer "{params.r_primer}" \\
            --locus "{params.locus}" \\
            --min-length {params.min_length} \\
            2>&1 | tee {log}

        # Validation: Check if output files were created and are not empty
        if [ ! -s "{output.fasta}" ] || [ ! -s "{output.taxonomy}" ]; then
            echo "Error: Read extraction resulted in empty files. Check primers and input data." >&2
            exit 1
        fi
        
        # Report the number of extracted records
        record_count=$(grep -c "^>" {output.fasta})
        echo "Successfully extracted $record_count records."
        """

# ---------------------------------------------------------------------
# Rule 2.3: Clean and validate FASTA and taxonomy files
# ---------------------------------------------------------------------
rule bold_clean_and_validate:
    """
    Clean and validate FASTA and taxonomy files before QIIME2 import.
    
    This is a critical QC step that mirrors the validation checks in Kim's RMD.
    It ensures that the sequence and taxonomy files are perfectly synchronized,
    removes sequences with invalid characters or that are too short, and handles
    any duplicate IDs. This prevents many common errors during QIIME2 import.
    """
    input:
        # Takes the raw extracted FASTA and taxonomy files.
        fasta = f"{workdir}/bold/bold_trimmed_seqs.fasta",
        taxonomy = f"{workdir}/bold/bold_trimmed_taxonomy.tsv"
    output:
        # Produces cleaned and validated versions of the files.
        cleaned_fasta = f"{workdir}/bold/bold_cleaned_seqs.fasta",
        cleaned_taxonomy = f"{workdir}/bold/bold_cleaned_taxonomy.tsv"
    params:
        # Minimum sequence length for a record to be kept.
        min_length = bold_min_length
    log:
        f"{workdir}/logs/bold_cleaning.log"
    conda:
        "../../envs/bold-pipeline.yml"
    shell:
        """
        set -euo pipefail
        
        echo "Cleaning and validating FASTA and taxonomy files..."
        python workflow/scripts/clean_bold_fasta_and_taxonomy.py \\
            --input-fasta {input.fasta} \\
            --input-taxonomy {input.taxonomy} \\
            --output-fasta {output.cleaned_fasta} \\
            --output-taxonomy {output.cleaned_taxonomy} \\
            --min-length {params.min_length} \\
            2>&1 | tee {log}

        # Validation: Check that output files are not empty and record counts match
        if [ ! -s "{output.cleaned_fasta}" ]; then
            echo "Error: Cleaned FASTA file is empty after validation script." >&2
            exit 1
        fi
        fasta_count=$(grep -c "^>" {output.cleaned_fasta})
        tax_count=$(wc -l < {output.cleaned_taxonomy})
        echo "$fasta_count records remain after cleaning and validation."
        if [ "$fasta_count" -ne "$tax_count" ]; then
            echo "Error: Mismatch between cleaned FASTA ($fasta_count) and taxonomy ($tax_count) records." >&2
            exit 1
        fi
        """ 