# Snakemake workflow: `colora`

[![Snakemake](https://img.shields.io/badge/snakemake-≥6.3.0-brightgreen.svg)](https://snakemake.github.io)
[![GitHub actions status](https://github.com/LiaOb21/colora/workflows/Tests/badge.svg?branch=main)](https://github.com/LiaOb21/colora/actions?query=branch%3Amain+workflow%3ATests)


A Snakemake workflow for for genome assembly.

![Alt text](Colora-1.jpg)
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

- [x] Rule for Nanoplot
- [x] Rule for fastp 
- [x] Rules for arima pipeline - split in several rules
- [x] Rule for yahs
- [ ] Formatting and linting to be fixed according to snakemake requirements
- [ ] log files: some of them are empty because it's impossible to redirect stderr and stdout to the file
- [ ] add optional params to all the rules
- [ ] implement ncbi `FCS` (decontamination) as optional rule (orange path in the scheme above)
- [ ] implement `assemblyQC` - waiting for new Merqury release to make a new conda recipe (light green path above)
- [ ] add singularity and docker as option for environment management (n.b. ncbi FCS can be run only with singularity)
- [ ] packages versions: I put the versions I used as mandatory, but I'm not sure if it is a good idea
- [ ] Rule `purge_dups.smk` and `purge_dups_alt.smk`: redirecting outputs to the final directory doesn't looks nice + in the root directory at the end of the workflow there are some files that I'm not sure why they are there 
- [ ] integrate the snakemake report in the workflow


Notes:
- purge dups rule wasn't working due to issues with the path where the output files are written.Evaluate possibility to change directory within the rule. Purge_dups rule wasn't working also because the command to convert the hifiasm gfa to fasta not interpreted correctly by snakemake. SOLVED


- Arima pipeline - changes compared to the original pipeline:
   - creating conda environments with needed tools so no need to specify tools' path
   - Remove the PREFIX line and the option -p $PREFIX from the bwa command, it is not necessary and creates problems in the reading of files
  - add -M flag in bwa mem command - step 1.A and 1.B
  - pipeline split in several rules

