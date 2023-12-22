# Snakemake workflow: `colora`

[![Snakemake](https://img.shields.io/badge/snakemake-≥6.3.0-brightgreen.svg)](https://snakemake.github.io)
[![GitHub actions status](https://github.com/LiaOb21/colora/workflows/Tests/badge.svg?branch=main)](https://github.com/LiaOb21/colora/actions?query=branch%3Amain+workflow%3ATests)


A Snakemake workflow for for genome assembly.

Input files: hifi reads and hic reads.

## Usage

The usage of this workflow is described in the [Snakemake Workflow Catalog](https://snakemake.github.io/snakemake-workflow-catalog/?usage=<LiaOb21>%2F<colora>).

If you use this workflow in a paper, don't forget to give credits to the authors by citing the URL of this (original) <colora> sitory and its DOI (see above).

- place raw hifi reads in /resources/raw_hifi
- place oatk database of interest from github.oatkdb.repo in /resources/oatkDB
- install oatk and make sure it's in your path


How to run `colora`:
```
conda config --set channel_priority flexible # if you want to use conda/mamba otherwise genomescope won't work

snakemake --software-deployment-method conda --conda-frontend mamba --snakefile workflow/Snakefile --cores all

snakemake --cluster "sbatch --cpus-per-task 50 --mem 200G -t 3-00" --use-conda --conda-frontend mamba --snakefile workflow/Snakefile
```

Before executing the command, ensure you have appropriately changed your `config.yaml`

# TODO


* The workflow will occur in the snakemake-workflow-catalog once it has been made public. Then the link under "Usage" will point to the usage instructions if `<owner>` and `<repo>` were correctly set.