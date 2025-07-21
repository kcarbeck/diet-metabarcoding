# File: workflow/rules/2_diet_analysis/4_assign_taxonomy.smk
# Description: classify ASVs and build final taxa barplot
# lower-case comments for consistency

project   = config["project_name"]
workdir   = f"results/{project}"

# where is the trained classifier? allow override in yaml; else default per-project path
classifier_qza = config.get("classifier_path", f"references/{project}_bold_final_classifier.qza") # ← must point at final alignment-based classifier
metadata_tsv   = config["metadata_tsv"]

qiime_env = "../../envs/qiime2-2025.4.yml"

# ---------------------------------------------------------------------
# 1) classify rep sequences – depends on classifier
#    classifier itself may be built by rules in rescript.smk (included separately)
# ---------------------------------------------------------------------

rule classify_taxonomy:
    input:
        rep_seqs_qza = f"{workdir}/dada2/rep_seqs.qza"
    params:
        classifier = classifier_qza
    output:
        taxonomy_qza = f"{workdir}/taxonomy/taxonomy.qza",
        taxonomy_qzv = f"{workdir}/taxonomy/taxonomy.qzv"
    log:
        f"{workdir}/logs/taxonomy.log"
    conda:
        qiime_env
    shell:
        (
            "set -euo pipefail; "
            "mkdir -p {workdir}/taxonomy {workdir}/logs; "
            # check classifier exists
            "if [ ! -f {params.classifier} ]; then echo 'ERROR: Classifier file {params.classifier} not found. Please build it first by running snakemake build_final_classifier.' >&2; exit 1; fi; "
            # classify sequences
            "qiime feature-classifier classify-sklearn "
            " --i-classifier {params.classifier} "
            " --i-reads {input.rep_seqs_qza} "
            " --o-classification {output.taxonomy_qza} "
            " 2>&1 | tee {log}; "
            # tabulate classification for easy viewing
            "qiime metadata tabulate "
            " --m-input-file {output.taxonomy_qza} "
            " --o-visualization {output.taxonomy_qzv} "
            " 2>&1 | tee -a {log}"
            "echo 'taxonomy step complete. taxonomy at: {output.taxonomy_qza}' | mail -s '[diet-metabarcoding] Step complete: taxonomy' {config[email]}"
        )

# (optional) unfiltered barplot for comparison
# -----------------------------------------------------------------------------

rule taxa_barplot_unfiltered:
    input:
        table_qza    = f"{workdir}/dada2/table.qza",
        taxonomy_qza = f"{workdir}/taxonomy/taxonomy.qza",
        metadata_tsv = metadata_tsv
    output:
        qzv = f"{workdir}/visualization/taxa_barplot_unfiltered.qzv"
    log:
        f"{workdir}/logs/barplot_unfiltered.log"
    conda:
        qiime_env
    shell:
        (
            "set -euo pipefail; mkdir -p {workdir}/visualization {workdir}/logs; "
            "qiime taxa barplot "
            " --i-table {input.table_qza} "
            " --i-taxonomy {input.taxonomy_qza} "
            " --m-metadata-file {input.metadata_tsv} "
            " --o-visualization {output.qzv} "
            " 2>&1 | tee {log}"
        )

# ---------------------------------------------------------------------
# 2) final stacked barplot using filtered table (after filter.smk)
# ---------------------------------------------------------------------

rule taxa_barplot:
    input:
        table_qza    = f"{workdir}/filter/filtered_table.qza",
        taxonomy_qza = f"{workdir}/taxonomy/taxonomy.qza",
        metadata_tsv = metadata_tsv
    output:
        barplot_qzv = f"{workdir}/visualization/taxa_barplot.qzv"
    log:
        f"{workdir}/logs/barplot.log"
    conda:
        qiime_env
    shell:
        (
            "set -euo pipefail; "
            "mkdir -p {workdir}/visualization {workdir}/logs; "
            "qiime taxa barplot "
            " --i-table {input.table_qza} "
            " --i-taxonomy {input.taxonomy_qza} "
            " --m-metadata-file {input.metadata_tsv} "
            " --o-visualization {output.barplot_qzv} "
            " 2>&1 | tee {log}"
        ) 