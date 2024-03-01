# Snakemake workflow: `colora`

[![Snakemake](https://img.shields.io/badge/snakemake-≥6.3.0-brightgreen.svg)](https://snakemake.github.io)
[![GitHub actions status](https://github.com/LiaOb21/colora/workflows/Tests/badge.svg?branch=main)](https://github.com/LiaOb21/colora/actions?query=branch%3Amain+workflow%3ATests)
[![DOI](https://zenodo.org/badge/730752023.svg)](https://zenodo.org/doi/10.5281/zenodo.10728679)

A Snakemake workflow for for genome assembly.

Why colora? :snake: Colora means "snake" in Sardinian language :snake: 

![Colora_1](https://github.com/LiaOb21/colora/assets/96196229/0d7cce80-bef9-46de-9213-e7e9343aa168)


Input reads: hifi reads, optionally ONT, and hic reads.
Other inputs: oatk database, ncbi FCS database (optional), BUSCO database (to be implemented)

## Usage

The usage of this workflow is described in the [Snakemake Workflow Catalog](https://snakemake.github.io/snakemake-workflow-catalog/?usage=LiaOb21%2Fcolora).

If you use this workflow in a paper, don't forget to give credits to the authors by citing the URL of this (original) <colora> sitory and its DOI (see above).

- place raw hifi reads in resources/raw_hifi
- place oatk database of interest from github.oatkdb.repo in resources/oatkDB
- place raw hic reads in resources/raw_hic
- place ncbi database for FCS-GX in resources/gx_db (optional, this needs ~500GB of disk space and a large RAM)


How to run `colora`:
```
snakemake --software-deployment-method conda --snakefile workflow/Snakefile --cores all

snakemake --software-deployment-method conda --snakefile workflow/Snakefile --cores all --dry-run


#for the cluster:

snakemake --software-deployment-method conda --conda-frontend mamba --snakefile workflow/Snakefile --cores 100
```

Before executing the command, ensure you have appropriately changed your `config.yaml`

Test the pipeline:

- 1. Download test data
- 2. Download oatk DB

```
git clone https://github.com/c-zhou/OatkDB.git
cd colora/resources
mkdir oatkDB
cp path/to/where/you/cloned/OatkDB/v20230921/dikarya* oatkDB/
```


- 3. Download FCS-GX test database 

You can skip this step if you are not going to run the decontamination step with FCS-GX
```
mamba create -n ncbi_fcsgx ncbi-fcs-gx
mamba activate ncbi_fcsgx
cd colora/resources
mkdir gx_test_db
cd gx_test_db
sync_files.py get --mft https://ftp.ncbi.nlm.nih.gov/genomes/TOOLS/FCS/database/test-only/test-only.manifest --dir ./test-only
```

- 4. Run the test pipeline

```
snakemake --configfile config/config_test.yaml --software-deployment-method conda --snakefile workflow/Snakefile --cores 4
```

# TODO


* The workflow will occur in the snakemake-workflow-catalog once it has been made public. Then the link under "Usage" will point to the usage instructions if `<owner>` and `<repo>` were correctly set.

- [x] Rule for Nanoplot
- [x] Rule for fastp 
- [x] Rules for arima pipeline - split in several rules
- [x] Rule for yahs
- [x] integrate the snakemake report in the workflow: not necessary
- [x] input / output: hardcoded is okay
- [x] test dataset
- [x] test config file
- [x] test possibility to add ONT reads as optional param in hifiasm
- [ ] test possibility to add HiC reads as optional params in hifiasm: file names change in this case. Need more study. Probably this needs a separate rule.
- [ ] packages versions: create stable yaml files with conda export
- [ ] add singularity and docker as option for environment management
- [x] implement ncbi `FCS` (decontamination) as optional rule (orange path in the scheme above)
- [x] make purging steps optional 
- [x] slurm integration (profile)
- [ ] setting of resources for each rule
- [x] Rule `purge_dups.smk` and `purge_dups_alt.smk`: redirecting outputs 
- [ ] implement `assemblyQC` - waiting for new Merqury release to make a new conda recipe (light green path above)
- [x] Formatting and linting to be fixed according to snakemake requirements
- [ ] log files: some of them are empty because it's impossible to redirect stderr and stdout to the file


Notes:

- Arima pipeline - changes compared to the original pipeline:
   - creating conda environments with needed tools so no need to specify tools' path
   - Remove the PREFIX line and the option -p $PREFIX from the bwa command, it is not necessary and creates problems in the reading of files
  - add -M flag in bwa mem command - step 1.A and 1.B
  - pipeline split in several rules

