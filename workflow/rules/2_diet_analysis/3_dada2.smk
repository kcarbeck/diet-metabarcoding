# File: workflow/rules/2_diet_analysis/3_dada2.smk
# Description: denoise paired-end reads with qiime dada2

# ---------------------------------------------------------------------
# config handles
# ---------------------------------------------------------------------
project   = config["project_name"]
workdir   = f"results/{project}"

trim_left_f = config.get("trim_left_f", 0)
trim_left_r = config.get("trim_left_r", 0)
trunc_len_f = config.get("trunc_len_f", 0)  # 0 means "do not truncate" in qiime
trunc_len_r = config.get("trunc_len_r", 0)
min_overlap = config.get("min_overlap", 12) # default in qiime is 12

# ---------------------------------------------------------------------
# rule: dada2 denoise-paired
# ---------------------------------------------------------------------

rule dada2_denoise:
    input:
        trimmed_qza = f"{workdir}/trim/trimmed_pe.qza"
    output:
        rep_seqs_qza = f"{workdir}/dada2/rep_seqs.qza",
        table_qza    = f"{workdir}/dada2/table.qza",
        stats_qza    = f"{workdir}/dada2/denoise_stats.qza",
        rep_seqs_qzv = f"{workdir}/dada2/rep_seqs.qzv",
        table_qzv    = f"{workdir}/dada2/table.qzv",
        stats_qzv    = f"{workdir}/dada2/denoise_stats.qzv"
    log:
        f"{workdir}/logs/dada2.log"
    conda:
        qiime_env
    threads:
        max(1, snakemake.threads)  # ensure >=1
    shell:
        (
            "set -euo pipefail; "
            "mkdir -p {workdir}/dada2 {workdir}/logs; "
            # run dada2 denoise-paired
            "qiime dada2 denoise-paired "
            " --i-demultiplexed-seqs {input.trimmed_qza} "
            " --p-trim-left-f {trim_left_f} "
            " --p-trim-left-r {trim_left_r} "
            " --p-trunc-len-f {trunc_len_f} "
            " --p-trunc-len-r {trunc_len_r} "
            " --p-min-overlap {min_overlap} "
            f" --p-n-threads {snakemake.threads} "
            " --o-representative-sequences {output.rep_seqs_qza} "
            " --o-table {output.table_qza} "
            " --o-denoising-stats {output.stats_qza} "
            " 2>&1 | tee {log}; "
            # summarise outputs for easy viewing
            "qiime feature-table summarize "
            " --i-table {output.table_qza} "
            " --o-visualization {output.table_qzv} "
            " 2>&1 | tee -a {log}; "
            "qiime feature-table tabulate-seqs "
            " --i-data {output.rep_seqs_qza} "
            " --o-visualization {output.rep_seqs_qzv} "
            " 2>&1 | tee -a {log}; "
            "qiime metadata tabulate "
            " --m-input-file {output.stats_qza} "
            " --o-visualization {output.stats_qzv} "
            " 2>&1 | tee -a {log}; "
            "echo 'DADA2 step complete. Inspect summaries for table, rep-seqs, and stats in {workdir}/dada2/' | mail -s '[diet-metabarcoding] Step complete: DADA2' {config[email]}"
        ) 