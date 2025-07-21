# File: workflow/rules/1_build_classifier/1_gather_bold.smk
# Description: Rule for gathering BOLD data, replicating the first part of Kim's pipeline.

# This rule downloads sequence data from the Barcode of Life Data System (BOLD)
# in taxonomic chunks to avoid API timeouts. It replicates the data gathering
# strategy from the beginning of `Kims_classifier_pipeline.rmd`.

# It also includes a critical feature: if a path to a raw, combined BOLD data
# file is provided in the configuration (`bold_raw_path`), the download step is
# skipped, and the provided file is symlinked into the results directory. This
# allows for rerunning the pipeline with a static, pre-existing dataset.

# ---------------------------------------------------------------------
# config handles
# ---------------------------------------------------------------------
project = config["project_name"]
workdir_path = f"results/{project}"


rule bold_gather_data:
    """
    Gather BOLD data from the API or use a pre-existing data file.
    
    This rule checks the user's configuration. If `bold_raw_path` is specified
    and points to a valid file, it creates a symbolic link to that file.
    Otherwise, it executes the `gather_bold_data.py` script to download all
    the necessary data from BOLD, following the taxonomic grouping strategy
    from Kim's RMD.
    """
    output:
        # The final output of this step is a single TSV file containing all the raw data.
        raw_data = f"{workdir_path}/bold/raw_bold_data.tsv"
    params:
        # The script uses this directory for its intermediate downloaded files.
        work_dir = f"{workdir_path}/bold",
        # The locus (e.g., 'COI-5P') is passed to the download script.
        locus = locus
    log:
        f"{workdir_path}/logs/bold_gather.log"
    conda:
        "../../envs/bold-pipeline.yml" # This environment contains the 'bold-py' library.
    shell:
        """
        set -euo pipefail
        mkdir -p {params.work_dir} {workdir_path}/logs
        
        # If a path to raw BOLD data is provided in the config, symlink it.
        # Otherwise, run the script to download the data from BOLD.
        if [ -n "{bold_raw_path}" ] && [ -f "{bold_raw_path}" ]; then
            echo "Using provided BOLD data: {bold_raw_path}"
            # Use absolute path for symlink to avoid issues
            ln -sf "$(realpath {bold_raw_path})" {output.raw_data}
        else
            echo "No pre-existing data found. Gathering comprehensive BOLD data..."
            # The python script manages its own subdirectories ('raw', 'processed')
            # and combines the final data into one file.
            python workflow/scripts/gather_bold_data.py \\
                --output-dir {params.work_dir} \\
                --locus "{params.locus}" \\
                2>&1 | tee {log}
            
            # The script outputs the final file to 'all_bold_data.tsv' inside its
            # working directory. We move it to the canonical output path for this rule.
            mv {params.work_dir}/all_bold_data.tsv {output.raw_data}
        fi
        """ 