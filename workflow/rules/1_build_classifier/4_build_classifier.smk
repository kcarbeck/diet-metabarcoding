# File: workflow/rules/4_build_classifier.smk
# Description: Rules for training and evaluating the BOLD taxonomic classifier.

# This file contains the final steps for building the BOLD classifier. It
# includes the main training rule and an optional evaluation rule to assess
# the performance of the trained classifier. This corresponds to the final
# part of Kim's original pipeline.

# ---------------------------------------------------------------------
# Rule 4.1: Train the BOLD classifier
# ---------------------------------------------------------------------
rule bold_train_classifier:
    """
    Train a Naive Bayes classifier from the dereplicated BOLD data.
    
    This is the primary rule for the classifier training. It uses the prepared
    reference sequences and taxonomy to train a scikit-learn Naive Bayes
    classifier, which is saved as a QIIME2 artifact. This is the main output
    of the first part of the pipeline.
    """
    input:
        # Depends on the dereplicated sequence and taxonomy artifacts.
        seqs_qza = f"{workdir}/bold/bold_derep_seqs.qza",
        tax_qza = f"{workdir}/bold/bold_derep_taxonomy.qza"
    output:
        # The final classifier artifact, saved to the references directory.
        classifier = classifier_qza
    params:
        # An email address can be provided in the config for notifications.
        email = config.get("email", "")
    log:
        f"{workdir}/logs/bold_train_classifier.log"
    conda:
        qiime_env
    shell:
        """
        set -euo pipefail
        mkdir -p references
        
        echo "Training BOLD Naive Bayes classifier..."
        qiime feature-classifier fit-classifier-naive-bayes \\
            --i-reference-reads {input.seqs_qza} \\
            --i-reference-taxonomy {input.tax_qza} \\
            --o-classifier {output.classifier} \\
            2>&1 | tee {log}
        
        # A simple validation check to ensure the output file was created.
        if [ ! -s "{output.classifier}" ]; then
            echo "Error: Classifier file '{output.classifier}' was not created or is empty." >&2
            exit 1
        fi

        # Send a completion email if an address was provided.
        if [ -n "{params.email}" ]; then
            echo "BOLD classifier has been built successfully at: {output.classifier}" | mail -s "[diet-metabarcoding] BOLD Classifier (Part 1) Complete" "{params.email}"
        fi
        """

# ---------------------------------------------------------------------
# Rule 4.2: Evaluate classifier performance (optional)
# ---------------------------------------------------------------------
rule bold_evaluate_classifier:
    """
    Evaluate the performance of the trained BOLD classifier.
    
    This optional rule uses RESCRIPt's evaluation tool to perform k-fold
    cross-validation on the reference dataset. It produces visualizations
    that help assess the accuracy and performance of the classifier at
    different taxonomic levels. This is computationally intensive and should
    be run deliberately.
    """
    input:
        # It requires the same inputs as the training rule, plus the trained classifier.
        seqs_qza = f"{workdir}/bold/bold_derep_seqs.qza",
        tax_qza = f"{workdir}/bold/bold_derep_taxonomy.qza",
        classifier = classifier_qza
    output:
        # The output is a directory containing QIIME2 visualization artifacts.
        evaluation_viz = directory(f"{workdir}/bold/classifier_evaluation")
    log:
        f"{workdir}/logs/bold_evaluate_classifier.log"
    conda:
        qiime_env
    shell:
        """
        set -euo pipefail
        
        echo "Evaluating BOLD classifier performance..."
        qiime rescript evaluate-fit-classifier \\
            --i-sequences {input.seqs_qza} \\
            --i-taxonomy {input.tax_qza} \\
            --i-classifier {input.classifier} \\
            --p-reads-per-batch 6000 \\
            --p-n-jobs 4 \\
            --output-dir {output.evaluation_viz} \\
            2>&1 | tee {log}
        """ 