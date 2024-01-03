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
- place raw hic reads in /resources/raw_hic


How to run `colora`:
```
conda config --set channel_priority flexible # if you want to use conda/mamba otherwise genomescope won't work

snakemake --software-deployment-method conda --conda-frontend mamba --snakefile workflow/Snakefile --cores all

snakemake --software-deployment-method conda --conda-frontend mamba --snakefile workflow/Snakefile --cores all --dry-run


#for the cluster:

srun --cpus-per-task 100 --mem 200G -t 5-00 snakemake --software-deployment-method conda --conda-frontend mamba --snakefile workflow/Snakefile --cores 100
```

Before executing the command, ensure you have appropriately changed your `config.yaml`

# TODO


* The workflow will occur in the snakemake-workflow-catalog once it has been made public. Then the link under "Usage" will point to the usage instructions if `<owner>` and `<repo>` were correctly set.

1. Rule for Nanoplot    DONE
2. Rule for fastp DONE - TO TEST WITH DATA
3. Rule for arima pipeline -> plan: make a directory with the perl scripts ued by the pipeline; place the main script of the pipeline in the rule; set the parameters required from the pipeline from the config file.
4. Rule for yahs
5. formatting
6. log files
7. add optional params to all the rules


Notes:
purge dups rule wasn't working due to issues with the path where the output files are written.Evaluate possibility to change directory within the rule.

Purge_dups rule wasn't working also because the command to convert the hifiasm gfa to fasta not interpreted correctly by snakemake.
Now hifiasm and purge_dups rules have been fixed.

Evaluate how to perform the second purge_dups run.

Arima pipeline:
 - creating conda environments with needed tools so no need to specify tools' path
 - Remove the PREFIX line and the option -p $PREFIX from the bwa command, it is not necessary and creates problems in the reading of files
- add -M flag in bwa mem command - step 1.A and 1.B

considera di mettere tutto in shell. No
splitting arima in several rules. Possible problem with the basename files of hic reads. to check