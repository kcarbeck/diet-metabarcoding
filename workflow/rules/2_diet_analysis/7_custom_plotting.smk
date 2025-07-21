# File: workflow/rules/2_diet_analysis/7_custom_plotting.smk
# Description: optional downstream visualisations in R
# these are placeholders; customise as you decide on final figures
# lower-case comments per lab style

project   = config["project_name"]
workdir   = f"results/{project}"

qiime_env = "../../envs/qiime2-2025.4.yml"  # for qiime tools export
r_env     = "../../envs/phyloseq.yml"        # create if you want an R-only env

# ---------------------------------------------------------------------
# rule: export filtered table + taxonomy as biom/tsv for R
# ---------------------------------------------------------------------

rule export_to_r:
    input:
        table_qza    = f"{workdir}/filter/filtered_table.qza",
        taxonomy_qza = f"{workdir}/taxonomy/taxonomy.qza"
    output:
        biom = f"{workdir}/export/feature-table.biom",
        tax  = f"{workdir}/export/taxonomy.tsv"
    log:
        f"{workdir}/logs/export.log"
    conda:
        qiime_env
    shell:
        (
            "set -euo pipefail; mkdir -p {workdir}/export {workdir}/logs; "
            "qiime tools export --input-path {input.table_qza} --output-path {workdir}/export 2>&1 | tee {log}; "
            "mv {workdir}/export/feature-table.biom {output.biom}; "
            "qiime tools export --input-path {input.taxonomy_qza} --output-path {workdir}/export 2>&1 | tee -a {log}; "
            "mv {workdir}/export/taxonomy.tsv {output.tax}"
        )

# ---------------------------------------------------------------------
# rule: create relative abundance barplot (example) via R script
# ---------------------------------------------------------------------

rule relabund_plot:
    input:
        biom = f"{workdir}/export/feature-table.biom",
        tax  = f"{workdir}/export/taxonomy.tsv",
        metadata = config["metadata_tsv"]
    output:
        pdf = f"{workdir}/visualization/rel_abund.pdf"
    log:
        f"{workdir}/logs/relabund_plot.log"
    conda:
        r_env
    shell:
        (
            "set -euo pipefail; mkdir -p {workdir}/visualization {workdir}/logs; "
            "Rscript workflow/scripts/rel_abund_plot.R {input.biom} {input.tax} {input.metadata} {output.pdf} 2>&1 | tee {log}"
        )

# note: provide the R script in workflow/scripts/rel_abund_plot.R (example uses phyloseq) 