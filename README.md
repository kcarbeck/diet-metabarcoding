# QIIME2 and BOLD Classifier Pipeline for Diet Metabarcoding

This repository contains a comprehensive Snakemake workflow for characterizing diets from fecal samples using QIIME2. It is designed to be modular, reproducible, and configurable, automating the entire process from raw sequencing reads to final analysis.

The workflow is divided into two main, independent pipelines:

1.  **Part 1. BOLD Classifier Construction (`build_final_classifier`)**: pipeline that downloads, filters, aligns, and processes data from the Barcode of Life Data System (BOLD) to create a custom, amplicon-specific taxonomic classifier. This is a heavyweight but essential preliminary step that only needs to be run once per primer set. This is adapted from Kim Navarro-Valez's R scripts and [Devon O'rourke's tuitorial](https://forum.qiime2.org/t/building-a-coi-database-from-bold-references/16129).
2.  **Part 2. QIIME2 Diet Analysis (`run_diet_pipeline`)**: The core analysis pipeline that takes your raw sequencing reads, performs QC (trimming, denoising), assigns taxonomy using the pre-built BOLD classifier, filters the data, and generates key visualizations and diversity analyses.

## File Structure

```
diet-metabarcoding/
├── README.md                          # you are here :p
├── Snakefile                          # main workflow orchestrator
├── config/
│   └── orchards.yaml                  # project-specific settings 
│   └── <yourproject>.yaml             # copy above and edit! **(EDIT ME)**
├── envs/
│   ├── bold-pipeline.yml              # conda environment for BOLD data handling
│   └── qiime2-2025.4.yml              # main QIIME2 Conda environment
├── workflow/
│   ├── rules/
│   │   ├── 1_build_classifier/      # rules for the BOLD Classifier Pipeline
│   │   │   ├── 1_gather_bold.smk
│   │   │   ├── 2_process_bold.smk
│   │   │   ├── 3_qiime_process.smk
│   │   │   ├── 4_build_classifier.smk
│   │   │   ├── 5_align_and_trim.smk
│   │   │   └── 6_finalize_classifier.smk
│   │   └── 2_diet_analysis/         # rules for the QIIME2 Diet Analysis Pipeline
│   │       ├── 1_import_and_validate.smk
│   │       ├── 2_cutadapt.smk
│   │       ├── 3_dada2.smk
│   │       ├── 4_assign_taxonomy.smk
│   │       ├── 5_filter_asvs.smk
│   │       ├── 6_diversity_and_rarefaction.smk
│   │       └── 7_custom_plotting.smk
│   └── scripts/                       # helper scripts used by the rules
│       ├── clean_bold_fasta_and_taxonomy.py
│       ├── extract_alignment_region.py
│       ├── extract_reads_bold.py
│       ├── filter_by_geo.py
│       ├── gather_bold_data.py
│       └── validate_metadata.py
├── data/
│   ├── raw/<project>/                 # gzipped FASTQ files (e.g., data/raw/orchards/)
│   └── metadata/                      # manifest and metadata files
├── references/                        # output directory for trained classifiers
└── results/                           # main output directory (auto-generated)
```

---

## How to Use

### 1. Initial Setup

**Clone the repository:**
```bash
git clone https://github.com/kcarbeck/diet-metabarcoding.git
cd diet-metabarcoding
```

**Create Conda environments:**
The workflow relies on two separate Conda environments. Create them using the provided YAML files. Using `mamba` is recommended for faster environment creation.
```bash
# It's recommended to use mamba for faster environment creation
conda install -n base -c conda-forge mamba
#otherwise can use 'conda' instead of 'mamba' below

#Download the specific QIIME2 environment file if you don't have it
# Note: The version should match the one specified in the rules.
# This example uses 2024.10, but adjust as needed.
# wget https://data.qiime2.org/distro/metagenome/qiime2-metagenome-2024.10-py310-linux-conda.yml \
#   -O envs/qiime2-metagenome-2024.10.yml

# create the environment for the BOLD pipeline
mamba env create -f envs/bold-pipeline.yml

# create the environment for the main QIIME2 analysis
mamba env create -f envs/qiime2-metagenome-2024.10.yml

# create snakemake env
mamba create -n snakemake snakemake -c conda-forge -c bioconda
mamba activate snakemake
```

### 2. Configure Your Project

Before running the pipeline, you must configure it for your project.

1.  **Copy the config template:**
    ```bash
    cp config/orchards.yaml config/my_project.yaml
    ```
2.  **Edit `config/my_project.yaml`:** Open the new file and carefully edit the parameters. Key settings include `project_name`, primer sequences, file paths, and filtering thresholds. See the comments in the file for detailed explanations of each parameter.

### 3. Build the BOLD Classifier (One-Time Step)

This is a one-time, computationally intensive job that builds the final, alignment-based classifier required for the diet analysis. It only needs to be run once per primer set.

The command below will execute the entire classifier construction pipeline:
```bash
snakemake --use-conda --configfile config/my_project.yaml --cores 8 build_final_classifier
```
This single command triggers a two-stage process (detailed in the Workflow Overview below):
1.  **Initial Data Processing:** Downloads, filters, and processes data from BOLD, training a preliminary classifier.
2.  **Alignment & Refinement:** Aligns the sequences to your specific primers, trims them precisely to the amplicon region, and trains the final, highly accurate classifier.

The final output will be saved to the path specified by `classifier_path` in your config file (e.g., `references/my_project_bold_final_classifier.qza`).

**Using Pre-Downloaded BOLD Data**
If you already have a comprehensive BOLD data file, you can skip the slow download step.
1.  Edit `config/my_project.yaml` and set the `bold_raw_path` parameter to the absolute path of your file.
2.  Run the `build_final_classifier` command as usual. The pipeline will automatically use your local file.

### 4. Prepare and Run the Diet Analysis Pipeline

Once the classifier exists, you can run the diet analysis on your samples.

**Prepare Data:**
1.  **Place raw FASTQ files** in a new directory inside `data/raw/`, e.g., `data/raw/my_project/`.
2.  **Create a QIIME2 manifest file** (`data/metadata/my_project_manifest.tsv`). See the "Helper Commands" section below for a command to generate this.
3.  **Create a sample metadata file** (`data/metadata/my_project_meta.tsv`).

**Run the Pipeline:**
```bash
snakemake --use-conda --configfile config/my_project.yaml --cores 8 run_diet_pipeline
```
This will run all steps from importing your raw reads to generating the final filtered and rarefied results.

---
## Workflow Overview & Rule Cheatsheet

### Part 1: BOLD Classifier Construction (`build_final_classifier`)

This pipeline builds a high-quality, alignment-based classifier from BOLD data, tailored to your specific primers. It proceeds in two main stages.

#### Stage A: Initial Data Processing & Reference Creation
This stage downloads and processes all the raw data from BOLD, culminating in a dereplicated set of sequences and taxonomies.

| Rule | Purpose |
| :--- | :--- |
| `bold_gather_data` | Downloads and merges all BOLD records by taxon, or links to a user-provided file. |
| `bold_geo_filter` | (Optional) Filters BOLD records by country, state, or bounding box. |
| `bold_extract_reads` | Extracts amplicon regions using your primers with in-silico PCR. |
| `bold_clean_and_validate` | QCs and synchronizes FASTA and taxonomy files, removing invalid/short sequences. |
| `bold_import_to_qiime` | Converts the cleaned text files to QIIME2 `.qza` artifacts. |
| `bold_dereplicate` | Dereplicates sequences and taxonomy using RESCRIPt to create a non-redundant reference set. |

#### Stage B: Alignment-Based Refinement
This stage takes the processed reference data and refines it to be specific to your amplicon region, creating a more accurate classifier. This is analogous to "Part 3" of the original R-based pipeline.

| Rule | Purpose |
| :--- | :--- |
| `bold_prepare_small_alignment_subsample` | Creates a high-quality subsample for building a reference alignment. |
| `bold_small_alignment` | Performs a small-scale MAFFT alignment on the subsample. |
| `bold_identify_primer_coordinates` | Maps primers to the small alignment to find exact trim coordinates. |
| `bold_prepare_large_alignment_set` | Prepares the full dataset by removing the sequences used in the small alignment. |
| `bold_large_alignment` | Aligns the full dataset to the small reference alignment using `mafft --addfull`. |
| `bold_extract_large_amplicon` | Trims the entire large alignment to the precise amplicon region identified earlier. |
| `bold_filter_and_import_final` | Imports, degaps, and length-filters the final amplicon sequences. |
| `bold_prepare_final_taxonomy` | Creates and validates the corresponding taxonomy for the final set of sequences. |
| `bold_final_dereplicate_and_train` | Performs a final dereplication and trains the final, alignment-based classifier. |

*For a detailed dependency graph, run `snakemake --dag | dot -Tpng > dag.png`.*


### Part 2: Diet Analysis Pipeline (`run_diet_pipeline`)

Once the classifier is built, this is the core pipeline for processing your metabarcoding data.

| Rule | Description | Key Output(s) |
| :--- | :--- | :--- |
| `run_diet_pipeline` | **(Target)** Runs the entire diet analysis from reads to plots. | `results/<project>/visualization/` |
| `validate_metadata` | Checks that manifest and metadata files are valid and consistent. | `logs/validate_metadata.done` |
| `qiime_import_demux` | Imports raw FASTQs and summarizes quality. | `demux/demux_pe.qzv` |
| `cutadapt_trim` | Trims primers and adapters from reads. | `trim/trimmed_pe.qzv` |
| `dada2_denoise` | Denoises reads into ASVs (Amplicon Sequence Variants). | `dada2/table.qzv` |
| `classify_taxonomy` | Assigns taxonomy to ASVs using the BOLD classifier. | `taxonomy/taxonomy.qzv` |
| `filter_table` | Filters the ASV table by taxonomy, abundance, and prevalence. | `filter/filtered_table.qzv` |
| `alpha_rarefaction` | Creates alpha rarefaction curves to guide sampling depth selection. | `diversity/alpha_rarefaction.qzv` |
| `rarefy_table` | Rarefies the feature table to a specified sampling depth. | `diversity/rarefied_table.qzv` |
| `taxa_barplot_rarefied` | Generates a barplot of the final, rarefied taxonomic composition. | `visualization/taxa_barplot_rarefied.qzv` |

---
<br>

## Helper Commands

**Useful Snakemake Arguments:**
*   `--dry-run` or `-n`: Show what would run without executing.
*   `--cores <N>`: Number of CPU cores to use.
*   `--use-conda`: Use the conda environments defined in the rule files.
*   `-p`: Print shell commands as they are executed.
*   `--rerun-incomplete`: Resume a failed or interrupted run.
*   `--until <rule>`: Run up to and including a specific rule.
*   `--forcerun <rule>`: Force a specific rule (and its children) to re-run.
*   `--keep-going`: Continue with independent jobs even if one fails.
<br>

**Inspecting Intermediate Steps**
You can run the pipeline up to a specific rule to inspect its output. For example, to stop after DADA2 to check the denoising stats:
```bash
# Run up to and including the dada2_denoise rule
snakemake --use-conda --configfile config/my_project.yaml --cores 8 --until dada2_denoise

# After inspecting the .qzv files in results/my_project/dada2/, you might update
# the trim/trunc parameters in your config file. Then, you can resume the
# pipeline, forcing DADA2 and all downstream steps to re-run:
snakemake --use-conda --configfile config/my_project.yaml --cores 8 --forcerun dada2_denoise run_diet_pipeline
``` 
<br>

**Generating the QIIME2 Manifest File:**

1. **Place raw FASTQ files in a new directory:** 
    `data/raw/<project>/` (e.g., `data/raw/orchards/`)
    <br>

2. **Create a QIIME2 manifest file** (`data/metadata/<project>_manifest.tsv`). This is a tab-separated file telling QIIME where to find your paired-end reads. **_Absolute paths are required!_**
    ```txt
    # Example: data/metadata/orchards_manifest.tsv

    sample-id   forward-absolute-filepath	                                                    reverse-absolute-filepath
    OR01        /path/to/your/workdir/diet-metabarcoding/data/raw/orchards/OR01_R1.fastq.gz     /path/to/your/workdir/diet-metabarcoding/data/raw/orchards/OR01_R2.fastq.gz
    OR02        /path/to/your/workdir/diet-metabarcoding/data/raw/orchards/OR02_R1.fastq.gz     /path/to/your/workdir/diet-metabarcoding/data/raw/orchards/OR02_R2.fastq.gz
    ...
    ```
    
    To help generate the required manifest file, run this command from the project's root directory. Remember to replace `<project>` with your actual project name (e.g., `orchards`).
    
    ```bash
    (cd data/raw/<project> && printf "sample-id\tforward-absolute-filepath\treverse-absolute-filepath\n" && for f in *_R1.fastq.gz; do id=${f%%_R1.fastq.gz}; printf "%s\t%s\t%s\n" "$id" "$(readlink -f "$f")" "$(readlink -f "${id}_R2.fastq.gz")"; done) > data/metadata/<project>_manifest.tsv
    ```

3. **Create a sample metadata file** (`data/metadata/<project>_meta.tsv`). This file contains information about your samples for downstream analysis. The pipeline will validate that the sample IDs in this file match the manifest.
   ```txt
   # Example: data/metadata/orchards_meta.tsv

   #SampleID   Host_Species    Plate           Type
   #q2:types   categorical     categorical     categorical
   OR01        SOSP            Plate1	    Sample
   OR_BLANK    none            Plate1          Blank
   ...
   ```

________________________

<br>

## Example workflow

### 1. Load snakemake

```bash
mamba activate snakemake
```

### 2. Update YAML
Open `config/<project>.yaml` and set parameters

If you have already downloaded BOLD database for training the classifier, be sure to set:
`bold_raw_path: "/absolute/path/to/all_bold_data.tsv"`

Adjust primers, locus, cores, email etc in the same file.

### 3. Running the pipeline (rule-by-rule)

```bash
# see every rule in the pipeline:
snakemake --list

# see every rule in classifier pipeline:
snakemake --list | grep '^bold_'

#draw a dependency graph for the classifier branch:
snakemake --dag build_final_classifier \
        | dot -Tpdf > bold_classifier_dag.pdf
```

Snakemake lets you stop after a given rule by using `--until` so you can stop and inspect results before moving on to the next rule.

Below are the command series for the BOLD classifier pipeline (assuming path pointing to previously downloaded BOLD data). Update number of cores!

```bash
############################################################################
# PART 1: quick “general” classifier (through dereplication and training)
############################################################################
snakemake --use-conda --cores 4 --until bold_gather_data            # raw_bold_data.tsv (symlink)
snakemake --use-conda --cores 4 --until bold_geo_filter             # bold_filtered.tsv
snakemake --use-conda --cores 4 --until bold_extract_reads          # bold_trimmed_seqs.fasta / taxonomy.tsv
snakemake --use-conda --cores 4 --until bold_clean_and_validate     # bold_cleaned_*.*
snakemake --use-conda --cores 4 --until bold_import_to_qiime        # bold_seqs.qza / taxonomy.qza
snakemake --use-conda --cores 4 --until bold_dereplicate            # bold_derep_*.qza
snakemake --use-conda --cores 4 --until bold_train_classifier       # references/<project>_bold_initial_classifier.qza


############################################################################
# PART 2: alignment workflow that trims to the exact amplicon
############################################################################
snakemake --use-conda --cores 4 --until bold_prepare_small_alignment_subsample
snakemake --use-conda --cores 4 --until bold_create_primers_fasta
snakemake --use-conda --cores 4 --until bold_small_alignment
snakemake --use-conda --cores 4 --until bold_identify_primer_coordinates
snakemake --use-conda --cores 4 --until bold_prepare_large_alignment_set
snakemake --use-conda --cores 4 --until bold_large_alignment
snakemake --use-conda --cores 4 --until bold_extract_large_amplicon
snakemake --use-conda --cores 4 --until bold_filter_and_import_final
snakemake --use-conda --cores 4 --until bold_prepare_final_taxonomy
snakemake --use-conda --cores 4 --until bold_final_dereplicate_and_train   #final classifier

# optional: evaluate the classifier? WORK IN PROGRESS
snakemake --use-conda --cores 4 --until bold_evaluate_classifier
```

**OTHER TIPS:**
- Use `snakemake -n --use-conda --until RULE` (“dry run”) first
  - this will print the shell commands without executing them
- All logs are written to `results/<project>/logs/`
  - open them in using nano or tail them live:
    ```bash
    tail -f results/<project>/logs/bold_geo_filter.log
    ```
- if you need to re-run a rule after changing parameters:
  ```bash 
   -R <rule_name> 
   #e.g.:
   snakemake -R bold_extract_reads 
   ```
   - This will force rerunning even if outputs already exist.
  
- Generating a DAG after each step (see above) is a great way to visualise progress and dependencies!