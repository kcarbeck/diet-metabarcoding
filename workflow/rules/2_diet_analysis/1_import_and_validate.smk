# File: workflow/rules/2_diet_analysis/1_import_and_validate.smk
# Description: manifest import + demux summary

# ---------------------------------------------------------------------
# config handles
# ---------------------------------------------------------------------
manifest     = config["samplesheet"]          # tsv with sample‑id + fastq paths
metadata_tsv = config.get("metadata_tsv", None)  # may be needed downstream
project      = config["project_name"]
workdir      = f"results/{project}"


# ---------------------------------------------------------------------
# rule: qiime import paired-end fastqs -> demux.qza + qzv summary
# ---------------------------------------------------------------------

rule validate_metadata:
    input:
        manifest,
        metadata_tsv 
    output:
        touch("results/{project}/logs/validate_metadata.done")
    log:
        f"{workdir}/logs/validate_metadata.log"
    conda:
        qiime_env
    shell:
        (
            "set -euo pipefail; mkdir -p {workdir}/logs; "
            "python workflow/scripts/validate_metadata.py {input.manifest} {input.metadata_tsv} 2>&1 | tee {log}; "
            "touch {output}"
        )

# update qiime_import_demux to require validate_metadata
rule qiime_import_demux:
    input:
        manifest,
        validate="results/{project}/logs/validate_metadata.done"
    output:
        demux_qza = f"{workdir}/demux/demux_pe.qza",
        demux_qzv = f"{workdir}/demux/demux_pe.qzv"
    params:
        # qiime expects absolute paths inside manifest – trust the user
        import_type   = "'SampleData[PairedEndSequencesWithQuality]'",
        import_format = "PairedEndFastqManifestPhred33V2"
    conda:
        qiime_env
    log:
        f"{workdir}/logs/qiime_import.log"
    shell:
        # break lines with \ to keep style neat
        (
            "set -euo pipefail; "
            "mkdir -p {workdir}/demux {workdir}/logs; "
            # qiime import
            "qiime tools import "
            " --type {params.import_type} "
            " --input-path {input} "
            " --input-format {params.import_format} "
            " --output-path {output.demux_qza} "
            " 2>&1 | tee {log}; "
            # qiime demux summarize
            "qiime demux summarize "
            " --i-data {output.demux_qza} "
            " --o-visualization {output.demux_qzv} "
            " 2>&1 | tee -a {log}; "
            "echo 'Import step complete. Demux summary at: {output.demux_qzv}' | mail -s '[diet-metabarcoding] Step complete: import' {config[email]}"
        ) 