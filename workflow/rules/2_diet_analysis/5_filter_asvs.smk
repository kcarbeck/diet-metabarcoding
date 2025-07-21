# File: workflow/rules/2_diet_analysis/5_filter_asvs.smk
# Description: drop unwanted features / taxa after dada2 + taxonomy

# ---------------------------------------------------------------------
# config handles
# ---------------------------------------------------------------------
project   = config["project_name"]
workdir   = f"results/{project}"

tax_include = config.get("tax_filter", "Arthropoda")  # regex include pattern
min_freq    = int(config.get("min_read_count", 0))
min_prev    = float(config.get("min_prevalence", 0))    # proportion (0–1)

# compute min number of samples from prevalence proportion
# re-parse manifest quickly to count samples (lightweight)
with open(config["samplesheet"], "r") as mf:
    _sample_n = len(mf.readlines()) - 1  # minus header
min_samples = int(round(min_prev * _sample_n)) if min_prev > 0 else 0

# helper to add param only if value > 0 or non-empty
_param = lambda flag, val: f" {flag} {val}" if val else ""

# ---------------------------------------------------------------------
# rule: filter feature table by taxonomy + freq/prevalence
# ---------------------------------------------------------------------

rule filter_table:
    input:
        table_qza   = f"{workdir}/dada2/table.qza",
        taxonomy_qza= f"{workdir}/taxonomy/taxonomy.qza"
    output:
        filt_table_qza = f"{workdir}/filter/filtered_table.qza",
        filt_table_qzv = f"{workdir}/filter/filtered_table.qzv"
    log:
        f"{workdir}/logs/filter.log"
    conda:
        qiime_env
    shell:
        (
            "set -euo pipefail; "
            "mkdir -p {workdir}/filter {workdir}/logs; "
            # 1) taxonomy include filter
            "qiime taxa filter-table "
            " --i-table {input.table_qza} "
            " --i-taxonomy {input.taxonomy_qza} "
            " --p-include '{tax_include}' "
            " --o-filtered-table _tmp_tax_filt.qza "
            " 2>&1 | tee {log}; "
            # 2) frequency / prevalence filter (skip flags if zero)
            "qiime feature-table filter-features "
            " --i-table _tmp_tax_filt.qza "
            "{_param('--p-min-frequency', min_freq)}"
            "{_param('--p-min-samples',  min_samples)}"
            " --o-filtered-table {output.filt_table_qza} "
            " 2>&1 | tee -a {log}; "
            # 3) summarise filtered table
            "qiime feature-table summarize "
            " --i-table {output.filt_table_qza} "
            " --m-sample-metadata-file {config[metadata_tsv]} "
            " --o-visualization {output.filt_table_qzv} "
            " 2>&1 | tee -a {log}; "
            "rm _tmp_tax_filt.qza; "
            "echo 'Feature table filtering complete. Inspect summary at: {output.filt_table_qzv}' | mail -s '[diet-metabarcoding] Step complete: Filtering' {config[email]}"
        ) 