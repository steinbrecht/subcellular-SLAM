# Subcellular SLAM seq: snakemake pipelines to fit subcellular T2C data
This repository contains three Snakemake pipelines to process metabolic labeling RNA sequencing data (fastq files) and fit a kinetic model to the processed data to estimate subcellular flow rates transcriptome-wide for each gene individually. To run the pipelines, you need to have conda installed. There are three main steps in the analysis pipeline:
- Genome alignment and intron genome generation
- T2C mutation counting and T2C data normalization
- Model fitting

### Genome alignment and intron genome generation
Corresponding folder: genome_alignment
This pipeline aligns read data to a genome and calculates transcript-per-million values. It also creates a genome annotation file that contains exclusively introns, which is used in the next step to calculate nuclear pre-mRNA data.

### T2C mutation counting and T2C data normalization
Corresponding folder: t2c_data_processing
This pipeline counts T2C mutations using the "CountT2C" program in the "t2c_count_program" folder. Then, it estimates conversion rates. Lastly, it estimates the ratio of new to total mRNA per intron and exon.

### Model fitting
Corresponding folder: model_fit
This pipeline fits a kinetic model to the fully-processed subcellular labeling data using the python package "lmfit". The are two models: 3-step and 4-step, each with its own Snakefile. Each Snakefile outputs two tables where each row corresponds to a specific gene. One table contains the best fit results for the kinetic parameters along with other fit statistics. The other table contains the varaiance of each kinetic parameter calculated from the ten best fit results.

## Corresponding raw data
The fastq files for which this pipeline was made can be found on GEO under accession number "GSE252199".

## Reproducibility
Instructions to run the the pipelines.
- install conda (or even better mamba) if necessary
- run `conda env create -f <placeholder_conda_env_file>` to create a new conda environment with all the necessary packages to run the code
- activate the environment with `conda activate <placeholder>`
- run Snakemake simply with `snakemake`
