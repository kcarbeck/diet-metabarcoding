# File: workflow/rules/1_build_classifier/3_qiime_process.smk
# Description: Rules for importing and processing data within the QIIME2 environment.

# This file handles the transition from standard text files (FASTA, TSV) into
# QIIME2 artifacts (.qza). It covers the import and dereplication steps, which
# are essential for preparing the reference data for classifier training.

# ---------------------------------------------------------------------
# config handles
# ---------------------------------------------------------------------
project = config["project_name"]
workdir_path = f"results/{project}"
forward_primer = config["forward_primer"]
reverse_primer = config["reverse_primer"]
locus = config["locus"]
bold_min_length = config["bold_min_length"]

# ---------------------------------------------------------------------
# Rule 3.1: Import sequences and taxonomy into QIIME2 format
# ---------------------------------------------------------------------
rule bold_import_to_qiime:
    """
    Import cleaned sequences and taxonomy into QIIME2 artifacts.
    
    This rule takes the validated FASTA and taxonomy files and converts them
    into the .qza format used by QIIME2. This is the gateway into the QIIME2
    portion of the pipeline.
    """
    input:
        # Depends on the cleaned and validated output from the previous step.
        fasta = f"{workdir_path}/bold/bold_cleaned_seqs.fasta",
        taxonomy = f"{workdir_path}/bold/bold_cleaned_taxonomy.tsv"
    output:
        # Produces two QIIME2 artifacts: one for sequences and one for taxonomy.
        seqs_qza = f"{workdir_path}/bold/bold_seqs.qza",
        tax_qza = f"{workdir_path}/bold/bold_taxonomy.qza"
    log:
        f"{workdir_path}/logs/bold_import.log"
    conda:
        "../../" + qiime_env # Uses the main QIIME2 conda environment.
    shell:
        """
        set -euo pipefail
        
        echo "Importing sequences to QIIME2..."
        qiime tools import \\
            --input-path {input.fasta} \\
            --output-path {output.seqs_qza} \\
            --type 'FeatureData[Sequence]' \\
            2>&1 | tee {log}
        
        echo "Importing taxonomy to QIIME2..."
        qiime tools import \\
            --input-path {input.taxonomy} \\
            --output-path {output.tax_qza} \\
            --type 'FeatureData[Taxonomy]' \\
            --input-format HeaderlessTSVTaxonomyFormat \\
            2>&1 | tee -a {log}
        """

# ---------------------------------------------------------------------
# Rule 3.2: Dereplicate sequences and taxonomy using RESCRIPt
# ---------------------------------------------------------------------
rule bold_dereplicate:
    """
    Dereplicate sequences and taxonomy to create a non-redundant reference set.
    
    This rule uses the `qiime rescript dereplicate` command to identify and
    collapse identical sequences. It intelligently handles taxonomy, ensuring that
    the resulting taxonomic annotations for the unique sequences are consistent.
    This step significantly reduces the size of the dataset and improves the
    quality of the final classifier.
    """
    input:
        # Takes the imported QIIME2 artifacts.
        seqs_qza = f"{workdir_path}/bold/bold_seqs.qza",
        tax_qza = f"{workdir_path}/bold/bold_taxonomy.qza"
    output:
        # Produces dereplicated versions of the sequence and taxonomy artifacts.
        derep_seqs = f"{workdir_path}/bold/bold_derep_seqs.qza",
        derep_tax = f"{workdir_path}/bold/bold_derep_taxonomy.qza"
    log:
        f"{workdir_path}/logs/bold_dereplicate.log"
    conda:
        "../../" + qiime_env
    shell:
        """
        set -euo pipefail
        
        echo "Dereplicating sequences and taxonomy with RESCRIPt..."
        qiime rescript dereplicate \\
            --i-sequences {input.seqs_qza} \\
            --i-taxa {input.tax_qza} \\
            --p-mode 'super' \\
            --p-derep-prefix \\
            --p-rank-handles kingdom phylum class order family genus species \\
            --o-dereplicated-sequences {output.derep_seqs} \\
            --o-dereplicated-taxa {output.derep_tax} \\
            2>&1 | tee {log}
        """ 