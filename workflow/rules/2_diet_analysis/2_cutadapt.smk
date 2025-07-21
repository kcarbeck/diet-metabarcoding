# File: workflow/rules/2_diet_analysis/2_cutadapt.smk
# Description: primer/adaptor trimming using qiime cutadapt trim‑paired

# ---------------------------------------------------------------------
# config handles (fall back to empty string if not provided)
# ---------------------------------------------------------------------

project   = config["project_name"]
workdir_path = f"results/{project}"

a_fwd_adapt = config.get("adapter_f", "")   # 3' adapter on fwd read
r_rev_adapt = config.get("adapter_r", "")   # 3' adapter on rev read
fwd_primer  = config.get("forward_primer", "")   # 5' primer on fwd read
rev_primer  = config.get("reverse_primer", "")   # 5' primer on rev read

# helper to build optional parameter strings; returns "" if seq is ""
def _param(flag, seq):
    return f" {flag} {seq}" if seq else ""

# ---------------------------------------------------------------------
# rule: cutadapt trim‑paired  (depends on qiime_import_demux output)
# ---------------------------------------------------------------------

rule cutadapt_trim:
    input:
        demux_qza = f"{workdir_path}/demux/demux_pe.qza"
    output:
        trimmed_qza = f"{workdir_path}/trim/trimmed_pe.qza",
        trimmed_qzv = f"{workdir_path}/trim/trimmed_pe.qzv"
    log:
        f"{workdir_path}/logs/cutadapt/cutadapt.log"
    conda:
        qiime_env
    shell:
        (
            "set -euo pipefail; "
            "mkdir -p {workdir_path}/trim {workdir_path}/logs; "
            # build cutadapt flags dynamically inside bash for readability
            "fwd_adapt='{_param('--p-adapter-f', a_fwd_adapt)}'; "
            "rev_adapt='{_param('--p-adapter-r', r_rev_adapt)}'; "
            "fwd_front='{_param('--p-front-f', fwd_primer)}'; "
            "rev_front='{_param('--p-front-r', rev_primer)}'; "
            # run cutadapt trim‑paired
            "qiime cutadapt trim-paired "
            " --i-demultiplexed-sequences {input.demux_qza} "
            "$fwd_adapt $rev_adapt $fwd_front $rev_front "
            " --p-discard-untrimmed False "
            " --o-trimmed-sequences {output.trimmed_qza} "
            " --verbose 2>&1 | tee {log}; "
            # summarize the trimmed sequences for quality check
            "qiime demux summarize "
            " --i-data {output.trimmed_qza} "
            " --o-visualization {output.trimmed_qzv} "
            " 2>&1 | tee -a {log}; "
            "echo 'cutadapt trim step complete. Trimmed sequence summary at: {output.trimmed_qzv}' | mail -s '[diet-metabarcoding] Step complete: cutadapt trim' {config[email]}"
        ) 