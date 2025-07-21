# File: workflow/rules/1_build_classifier/6_finalize_classifier.smk
# Description: Rules for finalizing the reference data and training the amplicon-specific classifier.

# This file contains the concluding steps of the alignment-based pipeline,
# mirroring the end of "Part 3" in Kim's RMD. It involves importing the
# trimmed sequences, filtering them by length, preparing the final corresponding
# taxonomy, performing a final dereplication, and training the ultimate
# amplicon-specific classifier.

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
# Rule 6.1: Filter by length and import final sequences to QIIME2
# ---------------------------------------------------------------------
rule bold_filter_and_import_final:
    """
    Filter final sequences by length, import to QIIME2, and degap.
    
    This rule takes the trimmed sequences from the large alignment, imports
    them into QIIME2, and then applies a minimum length filter after degapping.
    This is a critical step to remove sequences that are too short after
    trimming, ensuring only high-quality, full-length amplicons are used.
    """
    input:
        final_fasta = f"{workdir_path}/bold/alignment/ref_primers_trimmed_large.fasta"
    output:
        aligned_qza = f"{workdir_path}/bold/alignment/BOLD_anml_seqs_aligned.qza",
        degapped_qza = f"{workdir_path}/bold/alignment/BOLD_anml_seqs_aligned_nogaps.qza"
    params:
        min_length = config.get("bold_min_length", 165)
    log:
        f"{workdir_path}/logs/bold_filter_length_import.log"
    conda:
        "../../" + qiime_env
    shell:
        """
        set -euo pipefail
        
        echo "Importing final aligned sequences to QIIME2..."
        qiime tools import \\
            --input-path {input.final_fasta} \\
            --output-path {output.aligned_qza} \\
            --type 'FeatureData[AlignedSequence]' \\
            2>&1 | tee {log}
        
        echo "Degapping sequences and filtering by minimum length ({params.min_length} bp)..."
        qiime rescript degap-seqs \\
            --i-aligned-sequences {output.aligned_qza} \\
            --p-min-length {params.min_length} \\
            --o-degapped-sequences {output.degapped_qza} \\
            2>&1 | tee -a {log}
        
        echo "Length filtering and import complete."
        """

# ---------------------------------------------------------------------
# Rule 6.2: Prepare the final taxonomy file
# ---------------------------------------------------------------------
rule bold_prepare_final_taxonomy:
    """
    Prepare the final taxonomy file corresponding to the length-filtered sequences.
    """
    input:
        # We need the QIIME artifact of the sequences to get the final list of IDs
        degapped_qza = f"{workdir_path}/bold/alignment/BOLD_anml_seqs_aligned_nogaps.qza",
        # And the original full taxonomy file
        full_tax_tsv = f"{workdir_path}/bold/alignment/derep_taxonomy.tsv"
    output:
        final_taxonomy_qza = f"{workdir_path}/bold/alignment/ref_primers_large_taxonomy.qza"
    log:
        f"{workdir_path}/logs/bold_prepare_final_taxonomy.log"
    conda:
        "../../" + qiime_env
    shell:
        """
        set -euo pipefail
        tmp_dir="{workdir_path}/bold/alignment/tmp_final_tax"
        mkdir -p $tmp_dir

        echo "Exporting final sequence IDs to filter taxonomy..."
        qiime tools export --input-path {input.degapped_qza} --output-path $tmp_dir/seqs
        
        fasta_path="$tmp_dir/seqs/dna-sequences.fasta"
        grep '^>' $fasta_path | sed 's/^>//' > $tmp_dir/keep_ids.txt
        
        # Validation: Check sequence count
        seq_count=$(wc -l < $tmp_dir/keep_ids.txt)
        echo "Found $seq_count sequences to process for final taxonomy."
        if [ "$seq_count" -eq 0 ]; then
            echo "Error: No sequences remained after length filtering. Cannot proceed." >&2
            exit 1
        fi

        echo "Filtering and cleaning taxonomy..."
        # Filter the original taxonomy file to keep only the IDs present in our final sequence set
        awk 'FNR==NR{a[$1]; next} $1 in a' $tmp_dir/keep_ids.txt {input.full_tax_tsv} > $tmp_dir/filtered_taxonomy.tsv
        
        # Clean the taxonomy strings (e.g., remove extra whitespace and ensure tab separation)
        awk -F'\\t' 'BEGIN{{OFS="\\t"}} {{gsub(/[[:space:]]+/, "-", $2); print $1, $2}}' $tmp_dir/filtered_taxonomy.tsv > $tmp_dir/cleaned_taxonomy.tsv
        
        # Validation: Check that taxonomy count matches sequence count
        tax_count=$(wc -l < $tmp_dir/cleaned_taxonomy.tsv)
        echo "Found $tax_count corresponding taxonomy entries."
        if [ "$seq_count" -ne "$tax_count" ]; then
            echo "Error: Mismatch between final sequences ($seq_count) and taxonomy ($tax_count)." >&2
            exit 1
        fi
        
        echo "Importing final taxonomy to QIIME2..."
        qiime tools import \\
            --type 'FeatureData[Taxonomy]' \\
            --input-format HeaderlessTSVTaxonomyFormat \\
            --input-path $tmp_dir/cleaned_taxonomy.tsv \\
            --output-path {output.final_taxonomy_qza} \\
            2>&1 | tee {log}
        
        echo "Final taxonomy preparation complete."
        """

# ---------------------------------------------------------------------
# Rule 6.3: Final dereplication and classifier training
# ---------------------------------------------------------------------
rule bold_final_dereplicate_and_train:
    """
    Perform a final dereplication and train the amplicon-specific classifier.
    """
    input:
        degapped_qza = f"{workdir_path}/bold/alignment/BOLD_anml_seqs_aligned_nogaps.qza",
        taxonomy_qza = f"{workdir_path}/bold/alignment/ref_primers_large_taxonomy.qza"
    output:
        classifier = f"references/{project}_bold_final_classifier.qza"
    params:
        reads_per_batch = config.get("reads_per_batch", 6000),
        n_jobs = config.get("n_jobs", 6),
        email = config.get("email", "")
    log:
        f"{workdir_path}/logs/bold_train_final_classifier.log"
    conda:
        "../../" + qiime_env
    shell:
        """
        set -euo pipefail
        tmp_dir="{workdir_path}/bold/alignment/tmp_final_train"
        mkdir -p $tmp_dir

        echo "Performing final dereplication before training..."
        qiime rescript dereplicate \\
            --i-sequences {input.degapped_qza} \\
            --i-taxa {input.taxonomy_qza} \\
            --p-mode 'super' --p-derep-prefix \\
            --p-rank-handles kingdom phylum class order family genus species \\
            --o-dereplicated-sequences $tmp_dir/final_seqs.qza \\
            --o-dereplicated-taxa $tmp_dir/final_taxa.qza \\
            2>&1 | tee {log}

        echo "Training final amplicon-specific BOLD classifier..."
        qiime feature-classifier fit-classifier-naive-bayes \\
            --i-reference-reads $tmp_dir/final_seqs.qza \\
            --i-reference-taxonomy $tmp_dir/final_taxa.qza \\
            --o-classifier {output.classifier} \\
            2>&1 | tee -a {log}
        
        if [ ! -s "{output.classifier}" ]; then
            echo "Error: Final classifier file was not created or is empty." >&2
            exit 1
        fi
        
        if [ -n "{params.email}" ]; then
            echo "Final BOLD classifier (Part 2) has been trained: {output.classifier}" | mail -s "[diet-metabarcoding] BOLD Classifier (Part 2) Complete" "{params.email}"
        fi
        
        echo "Final classifier training complete."
        """ 