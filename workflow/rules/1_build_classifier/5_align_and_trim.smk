# File: workflow/rules/1_build_classifier/5_align_and_trim.smk
# Description: Rules for sequence alignment, coordinate identification, and amplicon trimming.

# This file encapsulates the logic from "Part 3" of `Kims_classifier_pipeline.rmd`.
# It details the process of creating a reference alignment from a small subsample,
# identifying primer coordinates, aligning the full dataset to the reference,
# and then trimming all sequences to the precise amplicon region. This ensures
# the final classifier is highly specific to the target primers.

# ---------------------------------------------------------------------
# Rule 5.1: Prepare subsample for small alignment
# ---------------------------------------------------------------------
rule bold_prepare_small_alignment_subsample:
    """
    Prepare a high-quality subsample for the initial small alignment.
    
    This rule exports the dereplicated sequences and taxonomy from the main
    pipeline, then applies a series of filters (length, taxonomy, and quality)
    to create a clean subset of insect sequences, mirroring Kim's approach.
    A random subsample of this set is then taken to create the reference for
    the initial alignment.
    """
    input:
        seqs_qza = f"{workdir}/bold/bold_derep_seqs.qza",
        tax_qza = f"{workdir}/bold/bold_derep_taxonomy.qza"
    output:
        subsample_fasta = f"{workdir}/bold/alignment/derep_subsample.fasta",
        full_seqs_fasta = f"{workdir}/bold/alignment/derep_seqs.fasta",
        full_tax_tsv = f"{workdir}/bold/alignment/derep_taxonomy.tsv"
    params:
        subsample_size = config.get("bold_subsample_size", 3000),
        min_len = 150,
        max_len = 1000
    log:
        f"{workdir}/logs/bold_prepare_subsample.log"
    conda:
        "../../envs/bold-pipeline.yml"
    shell:
        """
        set -euo pipefail
        tmp_dir="{workdir}/bold/alignment/tmp_subsample"
        mkdir -p $tmp_dir

        echo "Exporting QIIME2 artifacts for alignment..."
        qiime tools export --input-path {input.seqs_qza} --output-path $tmp_dir/seqs
        qiime tools export --input-path {input.tax_qza} --output-path $tmp_dir/tax
        
        mv $tmp_dir/seqs/dna-sequences.fasta {output.full_seqs_fasta}
        mv $tmp_dir/tax/taxonomy.tsv {output.full_tax_tsv}

        echo "Creating high-quality subset for alignment..."
        # Filter by length
        seqkit seq --min-len {params.min_len} --max-len {params.max_len} -w 0 {output.full_seqs_fasta} > $tmp_dir/len_filtered.fasta
        
        # Filter by taxonomy (insect species with full assignment)
        grep -i 'p__arthropoda.*c__insecta' {output.full_tax_tsv} | grep -v 's__$' | cut -f1 > $tmp_dir/insect_species.ids
        seqkit grep --pattern-file $tmp_dir/insect_species.ids $tmp_dir/len_filtered.fasta > $tmp_dir/taxa_filtered.fasta
        
        # Filter out ambiguous bases
        seqkit grep -s -i -p "[^acgt]" -v $tmp_dir/taxa_filtered.fasta > $tmp_dir/final_subset.fasta

        echo "Creating random subsample of {params.subsample_size} sequences..."
        seqkit sample --number {params.subsample_size} --rand-seed 101 $tmp_dir/final_subset.fasta > {output.subsample_fasta}
        
        echo "Subsample preparation complete."
        2>&1 | tee {log}
        """

# ---------------------------------------------------------------------
# Rule 5.2: Create primer sequences FASTA file
# ---------------------------------------------------------------------
rule bold_create_primers_fasta:
    """Create a FASTA file containing the forward and reverse primers."""
    output:
        primers_fasta = f"{workdir}/bold/alignment/anml_primers.fasta"
    params:
        f_primer = forward_primer,
        r_primer = reverse_primer
    log:
        f"{workdir}/logs/bold_create_primers.log"
    shell:
        """
        set -euo pipefail
        mkdir -p $(dirname {output.primers_fasta})
        
        echo "Creating primer FASTA file..."
        echo -e ">PRIMER_ANML-f\\n{params.f_primer}\\n>PRIMER_ANML-r\\n{params.r_primer}" > {output.primers_fasta}
        
        echo "Primer file created successfully."
        2>&1 | tee {log}
        """

# ---------------------------------------------------------------------
# Rule 5.3: Small alignment with MAFFT
# ---------------------------------------------------------------------
rule bold_small_alignment:
    """Perform a small alignment with MAFFT on the subsample."""
    input:
        subsample_fasta = f"{workdir}/bold/alignment/derep_subsample.fasta"
    output:
        aligned_fasta = f"{workdir}/bold/alignment/derep_subsample_aligned.fasta"
    params:
        threads = config.get("bold_mafft_threads", -1)
    log:
        f"{workdir}/logs/bold_small_alignment.log"
    conda:
        "../../envs/bold-pipeline.yml"
    shell:
        """
        set -euo pipefail
        echo "Starting small alignment with MAFFT..."
        mafft --thread {params.threads} --auto {input.subsample_fasta} > {output.aligned_fasta}
        echo "Small alignment complete."
        2>&1 | tee {log}
        """

# ---------------------------------------------------------------------
# Rule 5.4: Identify primer coordinates
# ---------------------------------------------------------------------
rule bold_identify_primer_coordinates:
    """Identify primer coordinates by aligning primers to the small alignment."""
    input:
        aligned_fasta = f"{workdir}/bold/alignment/derep_subsample_aligned.fasta",
        primers_fasta = f"{workdir}/bold/alignment/anml_primers.fasta"
    output:
        mapped_fasta = f"{workdir}/bold/alignment/ref_primers_aligned.fasta",
        map_file = f"{workdir}/bold/alignment/anml_primers.fasta.map"
    params:
        threads = config.get("bold_mafft_threads", -1)
    log:
        f"{workdir}/logs/bold_identify_coordinates.log"
    conda:
        "../../envs/bold-pipeline.yml"
    shell:
        """
        set -euo pipefail
        echo "Identifying primer coordinates with MAFFT..."
        mafft --multipair --addfragments {input.primers_fasta} --keeplength --thread {params.threads} --mapout --reorder {input.aligned_fasta} > {output.mapped_fasta}
        echo "Primer coordinate identification complete."
        2>&1 | tee {log}
        """

# ---------------------------------------------------------------------
# Rule 5.5: Prepare sequence set for large alignment
# ---------------------------------------------------------------------
rule bold_prepare_large_alignment_set:
    """
    Prepare the set of sequences for the large alignment.
    
    This rule takes the full set of dereplicated sequences and removes the
    sequences that were already used in the small reference alignment. This
    prevents re-aligning sequences that are already part of the reference.
    """
    input:
        full_seqs_fasta = f"{workdir}/bold/alignment/derep_seqs.fasta",
        subsample_fasta = f"{workdir}/bold/alignment/derep_subsample.fasta"
    output:
        large_alignment_set = f"{workdir}/bold/alignment/seqsForLargeAlignment.fasta"
    log:
        f"{workdir}/logs/bold_prepare_large_alignment_set.log"
    conda:
        "../../envs/bold-pipeline.yml"
    shell:
        """
        set -euo pipefail
        echo "Preparing sequences for large alignment by removing reference subsample..."
        
        # Extract IDs from the subsample to create a droplist
        grep '^>' {input.subsample_fasta} | sed 's/^>//' > {workdir}/bold/alignment/droplist.ids
        
        # Use seqkit to exclude the subsample sequences from the full set
        seqkit grep -v --pattern-file {workdir}/bold/alignment/droplist.ids {input.full_seqs_fasta} > {output.large_alignment_set}
        
        echo "Sequence set for large alignment is ready."
        2>&1 | tee {log}
        """

# ---------------------------------------------------------------------
# Rule 5.6: Large alignment with MAFFT
# ---------------------------------------------------------------------
rule bold_large_alignment:
    """Align the full dataset to the small reference alignment."""
    input:
        seqs_for_large_alignment = f"{workdir}/bold/alignment/seqsForLargeAlignment.fasta",
        ref_aligned = f"{workdir}/bold/alignment/ref_primers_aligned.fasta"
    output:
        large_aligned = f"{workdir}/bold/alignment/large_alignment.fasta"
    params:
        threads = config.get("bold_mafft_threads", -1)
    log:
        f"{workdir}/logs/bold_large_alignment.log"
    conda:
        "../../envs/bold-pipeline.yml"
    shell:
        """
        set -euo pipefail
        echo "Starting large alignment. This may take a significant amount of time and memory..."
        mafft --auto --addfull {input.seqs_for_large_alignment} --keeplength --thread {params.threads} {input.ref_aligned} > {output.large_aligned}
        echo "Large alignment complete."
        2>&1 | tee {log}
        """

# ---------------------------------------------------------------------
# Rule 5.7: Extract amplicon region from large alignment
# ---------------------------------------------------------------------
rule bold_extract_large_amplicon:
    """Extract the final amplicon region from the full alignment."""
    input:
        large_aligned = f"{workdir}/bold/alignment/large_alignment.fasta",
        map_file = f"{workdir}/bold/alignment/anml_primers.fasta.map"
    output:
        final_fasta = f"{workdir}/bold/alignment/ref_primers_trimmed_large.fasta"
    log:
        f"{workdir}/logs/bold_extract_large_amplicon.log"
    conda:
        "../../envs/bold-pipeline.yml"
    shell:
        """
        set -euo pipefail
        tmp_trimmed="{workdir}/bold/alignment/ref_anml_primers_trimmed.fasta"
        
        echo "Parsing trimming coordinates from MAFFT map file..."
        # Get the 1-based position of the last base of the forward primer
        start_pos=$(awk '/^>PRIMER_ANML-f/,/^>/{if ($0 ~ />PRIMER_ANML-f/ || $0 ~ /^#/) next; print}' {input.map_file} | tail -n 1 | awk '{{print $3}}')
        # Get the 1-based position of the first base of the reverse primer
        end_pos=$(awk '/^>PRIMER_ANML-r/,/^>/{if ($0 ~ />PRIMER_ANML-r/ || $0 ~ /^#/) next; print}' {input.map_file} | head -n 1 | awk '{{print $3}}')
        
        # The python script is 0-based, start-inclusive, end-exclusive.
        # To trim from 1-based (start_pos + 1) to (end_pos - 1), we calculate script params.
        py_start=$start_pos
        py_end=$((end_pos - 1))

        echo "Identified 1-based coordinates from map: Fwd primer ends at $start_pos, Rev primer starts at $end_pos"
        echo "Trimming to amplicon region using 0-based script coordinates: start=$py_start, end=$py_end"
        
        python workflow/scripts/extract_alignment_region.py \\
            -i {input.large_aligned} \\
            -o $tmp_trimmed \\
            -s $py_start \\
            -e $py_end
        
        echo "Removing primer headers and finalizing sequences..."
        grep -v -i ">PRIMER" $tmp_trimmed | seqkit seq --upper-case -w 0 > {output.final_fasta}

        # Final validation checks
        if grep -q -i ">PRIMER" {output.final_fasta}; then
            echo "Error: Primer sequences were found in the final trimmed FASTA file." >&2
            exit 1
        fi
        if [ ! -s "{output.final_fasta}" ]; then
            echo "Error: Final trimmed FASTA file is empty. Check alignment and coordinates." >&2
            exit 1
        fi

        echo "Large amplicon extraction and cleaning complete."
        2>&1 | tee -a {log}
        """ 