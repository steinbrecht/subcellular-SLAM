Instructions to run the model fit (assumes you are in a terminal in the directory of this file)

!!Note that only data for a couple of genes are included due to size reasons!!

1. Install snakemake in a conda environment. This can be done with the yaml file in the conda_envs folder via the command:
'conda env create --file conda-envs/snake-env.yaml'

2. activate the environment:
'conda activate snake'

3. run the Snakemake pipelines locally
'snakemake --snakefile Snakefile_4-step_model'
'snakemake --snakefile Snakefile_3-step_model'

Alternative command to run on your cluster
4. run the Snakemake pipelines on your cluster (this starts 900 parallel jobs)
'snakemake --snakefile Snakefile_4-step_model --profile brecht_profile/ -j 900 -R run_fit'
'snakemake --snakefile Snakefile_3-step_model --profile brecht_profile/ -j 900 -R run_fit'