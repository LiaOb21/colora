# Snakemake workflow: `colora`

[![Snakemake](https://img.shields.io/badge/snakemake-≥6.3.0-brightgreen.svg)](https://snakemake.github.io)
[![GitHub actions status](https://github.com/LiaOb21/colora/workflows/Tests/badge.svg?branch=main)](https://github.com/LiaOb21/colora/actions?query=branch%3Amain+workflow%3ATests)
[![DOI](https://zenodo.org/badge/730752023.svg)](https://zenodo.org/doi/10.5281/zenodo.10728679)

A Snakemake workflow for for genome assembly.

Why colora? :snake: Colora means "snake" in Sardinian language :snake: 

![Colora (1)](https://github.com/LiaOb21/colora/assets/96196229/483cf3e4-ceef-4b1a-846a-c0eec629189f)


Input reads: hifi reads, optionally ONT, and hic reads.
Other inputs: oatk database, ncbi FCS database (optional), BUSCO database (to be implemented)

## Usage

The usage of this workflow is described in the [Snakemake Workflow Catalog](https://snakemake.github.io/snakemake-workflow-catalog/?usage=LiaOb21%2Fcolora).

If you use this workflow in a paper, don't forget to give credits to the authors by citing the URL of this (original) <colora> sitory and its DOI (see above).

## Environment set up

You first have to install conda and snakemake.

## Reads

Reads must be in `fastq.gz` format.

If HiFi reads are in multiple files, these are automatically joined by colora. HiFi reads are not automatically filtered by colora, therefore if you want to filter them this must be done previously.

Hi-C reads are automatically filtered by colora using fastp. Fastp removes adapters for paired-end data with the parameter `--detect_adapter_for_pe`, which is always set in colora. If you are using Arima Hi-C library prep kit generated data, Arima mapping pipeline suggests to trim 5 bases from the 5' end of both read 1 and read 2, and this can be achieved automatically with colora, setting the right parameters in the config file (see [config/README.md](https://github.com/LiaOb21/colora/blob/main/config/README.md)).

ONT reads (if you have them) must be previously joined (if in multiple files) and filtered (if you want to).

## Other inputs

The other files that we need in order to run colora are:

- oatk database (mandatory)
- busco database (mandatory)
- ncbi FCS-gx database (optional, you can avoid this if you are not planning to automatically remove contaminants from your assembly)
- reference genome for the species under study (optional, you can use this if you want to compare your assembly with the reference genome with quast)

## Running colora

If everything is set up correctly and the `config.yaml` file has been updated according to your needs (see [config/README.md](https://github.com/LiaOb21/colora/blob/main/config/README.md)), you should be able to run `colora` with this simple command:

```
snakemake --software-deployment-method conda
```

N.B. if you are using a server where jobs are normally submitted through SLURM or other schedulers, you might consider setting up a snakemake profile in your system to handle job submission.


## Test the pipeline (update needed):

- 1. Download test data (they will be available soon)
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
- [x] test possibility to add HiC reads as optional params in hifiasm: file names change in this case. Need more study. Probably this needs a separate rule.
- [x] implement ncbi `FCS` (decontamination) as optional rule (orange path in the scheme above)
- [x] make purging steps optional 
- [x] slurm integration (profile)
- [x] setting of resources for each rule
- [x] Rule `purge_dups.smk` and `purge_dups_alt.smk`: redirecting outputs 
- [x] Formatting and linting to be fixed according to snakemake requirements
- [ ] implement `assemblyQC` - waiting for new Merqury release to make a new conda recipe (light green path above)
- [x] organelle QC
- [x] packages versions: create stable yaml files with conda export
- [ ] add singularity and docker as option for environment management
- [x] check if there is a better way to define oatk output
- [ ] set up GitHub actions



