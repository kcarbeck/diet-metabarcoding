# =================================================================
# Main Snakefile for the Diet Metabarcoding Pipeline
# =================================================================
# This file serves as the central orchestrator for the entire Snakemake workflow.
# It is responsible for loading the configuration, defining global variables,
# and including the modular rule files that constitute the pipeline's logic.

# -----------------------------------------------------------------
# Load Configuration and Define Global Variables
# -----------------------------------------------------------------
# Load the project-specific configuration file. This allows all
# parameters to be managed from a central YAML file.
configfile: "config/orchards.yaml"

# Define high-level variables that are used across multiple rules.
# This ensures consistency and makes the pipeline easier to maintain.
project = config["project_name"]
workdir: f"results/{project}"
qiime_env ="envs/qiime2-metagenome-2024.10.yml" # central QIIME environment

# Define variables used across multiple rule files
locus = config["locus"]
bold_raw_path = config.get("bold_raw_path", "")
classifier_qza = config.get("classifier_path", f"references/{project}_bold_final_classifier.qza")
workdir_path = f"results/{project}"

# BOLD classifier-specific parameters are now defined in the config
# and accessed directly within their respective rules. This keeps the
# main Snakefile cleaner and more focused on workflow orchestration.

# -----------------------------------------------------------------
# Include Modular Rule Files
# -----------------------------------------------------------------
# The pipeline is broken down into two main parts, each with its own set of rules.

# Part 1: Build a custom, alignment-based BOLD classifier.
include: "workflow/rules/1_build_classifier/1_gather_bold.smk"
include: "workflow/rules/1_build_classifier/2_process_bold.smk"
include: "workflow/rules/1_build_classifier/3_qiime_process.smk"
include: "workflow/rules/1_build_classifier/4_build_classifier.smk"
include: "workflow/rules/1_build_classifier/5_align_and_trim.smk"
include: "workflow/rules/1_build_classifier/6_finalize_classifier.smk"

# Part 2: Run the main QIIME2 diet analysis pipeline.
include: "workflow/rules/2_diet_analysis/1_import_and_validate.smk"
include: "workflow/rules/2_diet_analysis/2_cutadapt.smk"
include: "workflow/rules/2_diet_analysis/3_dada2.smk"
include: "workflow/rules/2_diet_analysis/4_assign_taxonomy.smk"
include: "workflow/rules/2_diet_analysis/5_filter_asvs.smk"
include: "workflow/rules/2_diet_analysis/6_diversity_and_rarefaction.smk"
include: "workflow/rules/2_diet_analysis/7_custom_plotting.smk"

# -----------------------------------------------------------------
# Global Handlers (e.g., for notifications)
# -----------------------------------------------------------------
onsuccess:
    """
    This block is executed upon successful completion of the entire pipeline.
    The main target rules also have their own completion notifications.
    """
    shell("echo 'Snakemake pipeline completed successfully.' | mail -s '[diet-metabarcoding] Pipeline Success' {config[email]}")

onerror:
    """
    This block is executed if any rule in the pipeline fails. It sends an
    email notification with the error details for long running jobs
    """
    shell("echo 'A step in the Snakemake pipeline failed. See logs for details.' | mail -s '[diet-metabarcoding] Pipeline FAILED' {config[email]}")


# -----------------------------------------------------------------
# Target Rules
# -----------------------------------------------------------------
# These are the main entry points for the user. They can be invoked
# from the command line, e.g., `snakemake build_final_classifier`.



rule build_initial_classifier:
    """
    Target to build the first, general-purpose BOLD classifier.
    """
    input:
        classifier_qza
    output:
        touch(f"{workdir_path}/bold/initial_classifier_complete.txt")
    params:
        email = config.get("email", "")
    shell:
        """
        echo "Initial BOLD classifier pipeline complete. Classifier at: {input}"
        if [ -n "{params.email}" ]; then
            echo "Initial BOLD classifier is complete. Find it at {input}" | mail -s "[diet-metabarcoding] Initial Classifier Complete" "{params.email}"
        fi
        touch {output}
        """

rule build_final_classifier:
    """
    Target to build the final, amplicon-specific classifier.
    """
    input:
        f"references/{project}_bold_final_classifier.qza"
    output:
        touch(f"{workdir_path}/bold/final_classifier_complete.txt")
    params:
        email = config.get("email", "")
    shell:
        """
        echo "Final, alignment-based classifier pipeline complete. Classifier at: {input}"
        if [ -n "{params.email}" ]; then
            echo "The final, alignment-based BOLD classifier is complete. Find it at {input}" | mail -s "[diet-metabarcoding] Final Classifier Complete" "{params.email}"
        fi
        touch {output}
        """

# -----------------------------------------------------------------
# Main Diet Analysis Pipeline Target
# -----------------------------------------------------------------

rule all:
    """
    Default target: runs the full diet analysis pipeline.
    This depends on the final classifier and produces all visualizations.
    """
    input:
        f"{workdir_path}/diet_analysis_complete.txt"


rule run_diet_pipeline:
    """
    Target to run the complete diet analysis from raw reads to final barplots.
    This is the main entry point for the diet analysis workflow.
    """
    input:
        f"{workdir_path}/visualization/taxa_barplot.qzv",
        f"{workdir_path}/visualization/taxa_barplot_unfiltered.qzv",
        f"{workdir_path}/visualization/taxa_barplot_rarefied.qzv",
        f"{workdir_path}/diversity/alpha_rarefaction.qzv",
        f"{workdir_path}/visualization/rel_abund.pdf",
        f"{workdir_path}/visualization/rarefied_rel_abund.pdf"
    output:
        touch(f"{workdir_path}/diet_analysis_complete.txt")
    params:
        email = config.get("email", "")
    shell:
        """
        echo "QIIME2 diet analysis pipeline is complete. Final visualizations are in {workdir_path}/visualization/ and {workdir_path}/diversity/"
        if [ -n "{params.email}" ]; then
            echo "The QIIME2 diet analysis pipeline is complete. Find final reports in {workdir_path}/visualization/ and {workdir_path}/diversity/" | mail -s "[diet-metabarcoding] Diet Analysis Complete" "{params.email}"
        fi
        """