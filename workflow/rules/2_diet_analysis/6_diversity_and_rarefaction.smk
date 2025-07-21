# File: workflow/rules/2_diet_analysis/6_diversity_and_rarefaction.smk
# Description: taxa collapse, alpha‑rarefaction, rarefy + extra barplots

project   = config["project_name"]
workdir   = f"results/{project}"
metadata_tsv = config["metadata_tsv"]
qiime_env = "../../envs/qiime2-2025.4.yml"

# yaml‑tunable params with legacy defaults
collapse_level        = int(config.get("collapse_level", 7))
alpha_min_depth       = int(config.get("alpha_min_depth", 500))
alpha_max_depth       = int(config.get("alpha_max_depth", 30000))
rarefy_sampling_depth = int(config.get("rarefy_sampling_depth", 2000))

# -----------------------------------------------------------------------------
# collapse filtered table to chosen taxonomic level
# -----------------------------------------------------------------------------

rule taxa_collapse:
    input:
        table_qza    = f"{workdir}/filter/filtered_table.qza",
        taxonomy_qza = f"{workdir}/taxonomy/taxonomy.qza"
    output:
        collapsed_qza = f"{workdir}/diversity/table_collapsed.qza"
    log:
        f"{workdir}/logs/taxa_collapse.log"
    conda:
        qiime_env
    shell:
        (
            "set -euo pipefail; mkdir -p {workdir}/diversity {workdir}/logs; "
            "qiime taxa collapse "
            " --i-table {input.table_qza} "
            " --i-taxonomy {input.taxonomy_qza} "
            " --p-level {collapse_level} "
            " --o-collapsed-table {output.collapsed_qza} "
            " 2>&1 | tee {log}"
        )

# -----------------------------------------------------------------------------
# alpha‑rarefaction plot on collapsed table
# -----------------------------------------------------------------------------

rule alpha_rarefaction:
    input:
        collapsed_qza = f"{workdir}/diversity/table_collapsed.qza",
        metadata_tsv  = metadata_tsv
    output:
        qzv = f"{workdir}/diversity/alpha_rarefaction.qzv"
    log:
        f"{workdir}/logs/alpha_rarefaction.log"
    conda:
        qiime_env
    shell:
        """
        qiime diversity alpha-rarefaction \
            --i-table {input.collapsed_qza} \
            --m-metadata-file {input.metadata_tsv} \
            --p-min-depth {alpha_min_depth} \
            --p-max-depth {alpha_max_depth} \
            --o-visualization {output.qzv} \
            --verbose &> {log}

        echo "Alpha rarefaction analysis complete. Visualization at: {output.qzv}" | mail -s "[diet-metabarcoding] Step complete: Diversity" {config[email]}
        """

# -----------------------------------------------------------------------------
# rarefy filtered table to fixed depth
# -----------------------------------------------------------------------------

rule rarefy_table:
    input:
        table_qza = f"{workdir}/filter/filtered_table.qza"
    output:
        rarefied_qza = f"{workdir}/diversity/rarefied_table.qza",
        rarefied_qzv = f"{workdir}/diversity/rarefied_table.qzv"
    log:
        f"{workdir}/logs/rarefy_table.log"
    conda:
        qiime_env
    shell:
        (
            "set -euo pipefail; mkdir -p {workdir}/diversity {workdir}/logs; "
            "qiime feature-table rarefy "
            " --i-table {input.table_qza} "
            " --p-sampling-depth {rarefy_sampling_depth} "
            " --o-rarefied-table {output.rarefied_qza} "
            " 2>&1 | tee {log}; "
            # summarise rarefied table
            "qiime feature-table summarize "
            " --i-table {output.rarefied_qza} "
            " --o-visualization {output.rarefied_qzv} "
            " 2>&1 | tee -a {log}"
        )

# -----------------------------------------------------------------------------
# barplot of rarefied table (final visual)
# -----------------------------------------------------------------------------

rule taxa_barplot_rarefied:
    input:
        table_qza    = f"{workdir}/diversity/rarefied_table.qza",
        taxonomy_qza = f"{workdir}/taxonomy/taxonomy.qza",
        metadata_tsv = metadata_tsv
    output:
        qzv = f"{workdir}/visualization/taxa_barplot_rarefied.qzv"
    log:
        f"{workdir}/logs/barplot_rarefied.log"
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