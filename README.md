# Subcellular SLAM seq: snakemake pipelines to fit subcellular T2C data
This repository contains three Snakemake pipelines to process metabolic labeling RNA sequencing data (fastq files) and fit a kinetic model to the processed data to estimate subcellular flow rates transcriptome-wide for each gene individually. This workflow was created while working on the paper [Subcellular mRNA kinetic modeling reveals nuclear retention as rate-limiting](https://www.biorxiv.org/content/10.1101/2024.03.11.584215v1) by Steinbrecht et al. There are three main steps in the workflow:
- Transcript quantification (TPM normalized)
- T2C mutation counting and T2C data normalization
- Model fitting

### Transcript quantification (TPM normalized)
*Corresponding folder*: `tpm_quantification`\
This pipeline aligns read data to a genome and calculates transcript-per-million values using `RSEM`. TPM data from all samples are combined in a single file `tpm_normalized_fraction_counts.csv`. The mouse genome (GENCODE M14, same used in paper) is autmatically downloaded and used as reference geome. Configure the pipeline to your need in the `config.yaml` file.

### T2C mutation counting and T2C data normalization
*Corresponding folder*: `slam_seq_analysis`\
This pipeline contains a full-scale SLAM-seq analysis (in a *GRAND-SLAM* manner). First, sequencing data is aligned using `STAR`. Then, labeling-induced T2C mutations are counted using the custom-built `CountT2C` program, where *T2C* SNPs are excluded. Next, conversion rates are estimated based on the distribution of *T2C* mutation frequencies over the number of *T* counts in a read. Lastly, the ratio of new to total mRNA per intron and exon is estimated by normalizing the *T2C* counts with the conversion rates. The mouse genome (GENCODE M14, same used in paper) is autmatically downloaded and used as reference geome. Configure the pipeline to your need in the `config.yaml` file.

### Model fitting
*Corresponding folder*: `model_fit`\
This pipeline fits a kinetic model to the fully-processed subcellular labeling data using the python package `lmfit`. The are two models: 3-step and 4-step, each with its own Snakefile. Each Snakefile outputs two tables where each row corresponds to a specific gene. One table contains the best fit results for the kinetic parameters along with other fit statistics. The other table contains the varaiance of each kinetic parameter calculated from the ten best fit results.

## Reproducibility
Instructions to run the the pipelines. The pipeline was tested on a SLURM-based high performance cluster.
- install conda (preferably [miniforge](https://github.com/conda-forge/miniforge))
- run `conda env create -n snake -c bioconda snakemake=7.*` to create a new conda environment with snakemake installed
- activate the environment with `conda activate snake`
- change the directory to the pipeline you want to run (e.g. `cd slam_seq_analysis` to do the SLAM-seq analysis)
- run the Snakemake pipleine with `snakemake`. The specific command to call snakemake can be found in the third line of each Snakefile. 

### Corresponding raw data
The fastq files for which this pipeline was made can be found on GEO under accession number [GSE252199](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE252199).

### Notes
The folder `brecht_profile` contains a SLURM profile that can be used on SLURM-based HP clusters to execute the snakemake workflows.
