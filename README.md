# Snakemake workflow: `colora`

[![Snakemake](https://img.shields.io/badge/snakemake-≥6.3.0-brightgreen.svg)](https://snakemake.github.io)
[![GitHub actions status](https://github.com/LiaOb21/colora/workflows/Tests/badge.svg?branch=main)](https://github.com/LiaOb21/colora/actions?query=branch%3Amain+workflow%3ATests)
[![DOI](https://zenodo.org/badge/730752023.svg)](https://zenodo.org/doi/10.5281/zenodo.10728679)

A Snakemake workflow for *de novo* genome assembly.

Why colora? :snake: Colora means "snake" in Sardinian language :snake: 

![Colora](https://github.com/LiaOb21/colora/assets/96196229/20fbc901-ae09-4c02-9278-f996e834a7ea)

## Overview

The aim of colora is to produce complete, chromosome-scale primary or phased assemblies by integrating the following tools:

- [Hifiasm](https://github.com/chhylp123/hifiasm): used to extract contigs from raw PacBio HiFi reads. Hifiasm operates either with HiFi reads exclusively or in conjunction with Oxford Nanopore reads to generate primary assemblies. Colora supports also the hifiasm 'Hi-C mode', used to create phased assemblies with distinct haplotypes.
- [FCS-GX](https://github.com/ncbi/fcs-gx): this pipeline is employed to eliminate contaminants from genome assemblies. This step is optional.
- [purge_dups](https://github.com/dfguan/purge_dups): applied to remove haplotypic duplications and overlaps from primary assemblies. This step is optional and must be skipped in case of phased assembly (Hi-C mode).
- [Arima Genomics Mapping Pipeline](https://github.com/ArimaGenomics/mapping_pipeline): used to map Hi-C reads to contigs. It has been adapted to Snakemake within Colora with minor modifications (`-M` flag added to `bwa mem` commands).
- [YaHS](https://github.com/c-zhou/yahs): used for scaffolding the assemblies.

In addition, Colora executes the following tasks:

- Quality assessment of raw HiFi reads using [NanoPlot](https://github.com/wdecoster/NanoPlot).
- Quality evaluation and filtering of Hi-C reads with [Fastp](https://github.com/OpenGene/fastp).
- Analysis of the k-mer spectrum from raw HiFi reads using [KMC](https://github.com/refresh-bio/KMC) and [GenomeScope2](https://github.com/tbenavi1/genomescope2.0).
- Assembly of mitochondria and chloroplasts (when applicable) with [Oatk](https://github.com/c-zhou/oatk).
- Quality evaluation of organelle assemblies with [gfastats](https://github.com/vgl-hub/gfastats) and [Bandage](https://github.com/rrwick/Bandage).
- Quality assessment of the *de novo* genome assembly throughout the workflow using [QUAST](https://github.com/ablab/quast) and [BUSCO](https://gitlab.com/ezlab/busco/-/releases#5.7.0).

## Usage

The usage of this workflow is described in the [Snakemake Workflow Catalog](https://snakemake.github.io/snakemake-workflow-catalog/?usage=LiaOb21%2Fcolora).

If you use this workflow in a paper, don't forget to give credits to the authors by citing the URL of this (original) <colora> sitory and its DOI (see above).

### 1. Environment set up

You first have to install conda/mamba and snakemake.

### 2. Clone colora git repo

```
git clone https://github.com/LiaOb21/colora.git
cd colora
mkdir resources # we link here all the input files in the following examples
```

### 3. Input reads

Reads must be in `fastq.gz` format.

If HiFi reads are in multiple files, these are automatically joined by colora. HiFi reads are not automatically filtered by colora, therefore if you want to filter them this must be done previously.

Hi-C reads are automatically filtered by colora using fastp. Fastp removes adapters for paired-end data with the parameter `--detect_adapter_for_pe`, which is always set in colora. If you are using Arima Hi-C library prep kit generated data, Arima mapping pipeline suggests to trim 5 bases from the 5' end of both read 1 and read 2, and this can be achieved automatically with colora, setting the right parameters in the config file (see [config/README.md](https://github.com/LiaOb21/colora/blob/main/config/README.md)).

ONT reads (if you have them) must be previously joined (if in multiple files) and filtered (if you want to).

### 4. Other inputs

The other files that we need in order to run colora are:

- **oatk database (mandatory)**

Dikarya in the following code is an example, you must choose the right database for your species:

```
git clone https://github.com/c-zhou/OatkDB.git
cd colora/resources
mkdir oatkDB
cd oatkDB
ln -s ~/software/OatkDB/v20230921/dikarya_mito.fam
ln -s ~/software/OatkDB/v20230921/dikarya_mito.fam.h3f
ln -s ~/software/OatkDB/v20230921/dikarya_mito.fam.h3i
ln -s ~/software/OatkDB/v20230921/dikarya_mito.fam.h3m
ln -s ~/software/OatkDB/v20230921/dikarya_mito.fam.h3p
```

- **busco database (mandatory)**

Go to https://busco-data.ezlab.org/v5/data/lineages/ and download the database of interest.
Even in this case, the following is just an example:

```
cd colora/resources
mkdir busco_db
cd busco_db
wget https://busco-data.ezlab.org/v5/data/lineages/saccharomycetes_odb10.2024-01-08.tar.gz
```

- **ncbi FCS-gx database (optional)** 

You can avoid this if you are not planning to automatically remove contaminants from your assembly.

```
mamba create -n ncbi_fcsgx ncbi-fcs-gx
mamba activate ncbi_fcsgx
cd colora/resources
mkdir gxdb
cd gxdb
sync_files.py get --mft https://ftp.ncbi.nlm.nih.gov/genomes/TOOLS/FCS/database/latest/all.manifest --dir ./gxdb
```
- **reference genome (optional)**

You can use the reference genome for the species under study if you want to compare your assembly with the reference genome with quast. You can also use the reference GFF file (see [config/README.md](https://github.com/LiaOb21/colora/blob/main/config/README.md)).

### 5. Running colora

If everything is set up correctly and the `config.yaml` file has been updated according to your needs (see [config/README.md](https://github.com/LiaOb21/colora/blob/main/config/README.md)), you should be able to run `colora` with this simple command:

```
snakemake --software-deployment-method conda
```

**Note:** if you are using a server where jobs are normally submitted through SLURM or other schedulers, you might consider setting up a snakemake profile in your system to handle job submission.


## Test the pipeline:

The test dataset provided is a subset of reads of the organism *Saccharomyces cerevisiae*. The data come from two different BioProjects:
- HiFi and ONT reads come from the BioProject [PRJNA1075684](https://www.ncbi.nlm.nih.gov/bioproject?LinkName=sra_bioproject&from_uid=31877222) (strain SPSC01)
- Hi-C reads come from the BioProject [PRJNA1013711](https://www.ncbi.nlm.nih.gov/bioproject?LinkName=sra_bioproject&from_uid=31374389) (strain YBP2)

This dataset is not supposed to have biological meaning, it has been crated only with the purpose of testing the workflow functionality. 

#### 1. Clone colora repository:

```
git clone https://github.com/LiaOb21/colora.git
cd colora
```

#### 2. Download oatk DB

```
git clone https://github.com/c-zhou/OatkDB.git
cd test_data
mkdir oatkDB
cd oatkDB
ln -s ~/software/OatkDB/v20230921/dikarya_mito.fam
ln -s ~/software/OatkDB/v20230921/dikarya_mito.fam.h3f
ln -s ~/software/OatkDB/v20230921/dikarya_mito.fam.h3i
ln -s ~/software/OatkDB/v20230921/dikarya_mito.fam.h3m
ln -s ~/software/OatkDB/v20230921/dikarya_mito.fam.h3p
cd ..
```
#### 3. Download busco lineage
  
```
mkdir busco_db
cd busco_db
wget https://busco-data.ezlab.org/v5/data/lineages/saccharomycetes_odb10.2024-01-08.tar.gz
tar -xzf saccharomycetes_odb10.2024-01-08.tar.gz
cd ..
```

#### 4. Download FCS-GX test database 

You can skip this step if you are not going to run the decontamination step with FCS-GX, in which case you should modify the `config/config_test.yaml` file setting `include_fcsgx: False`.

```
mamba create -n ncbi_fcsgx ncbi-fcs-gx
mamba activate ncbi_fcsgx
mkdir gx_test_db
cd gx_test_db
sync_files.py get --mft https://ftp.ncbi.nlm.nih.gov/genomes/TOOLS/FCS/database/test-only/test-only.manifest --dir ./test-only
mamba deactivate
cd ..
```

#### 5. Concatenate HiFi and ONT files

**Note:** With real data this step is only necessary for ONT reads (when available). Hifi files are automatically joined by the workflow. In this case, we have to perform this step manually because of the way the files are split. 


```
cd raw_hifi

cat hifi_test_SPSC01_SRR27947616_PRJNA1075684aa.fastq.gz hifi_test_SPSC01_SRR27947616_PRJNA1075684ab.fastq.gz > hifi_test_SPSC01_SRR27947616_PRJNA1075684.fastq.gz

rm hifi_test_SPSC01_SRR27947616_PRJNA1075684a*

cd ../raw_ont 

cat ont_test_SPSC01_SRR27947616_PRJNA1075684aa.fastq.gz ont_test_SPSC01_SRR27947616_PRJNA1075684ab.fastq.gz ont_test_SPSC01_SRR27947616_PRJNA1075684ac.fastq.gz > ont_test_SPSC01_SRR27947616_PRJNA1075684.fastq.gz

rm ont_test_SPSC01_SRR27947616_PRJNA1075684a*
cd ../..
```
 
#### 6. Run the test pipeline

```
mamba activate snakemake
snakemake --configfile config/config_test.yaml --software-deployment-method conda --snakefile workflow/Snakefile --cores 4
```

**Note:** The testing will take approximately 40 minutes. It may take longer depending on the time required for the downloading of the conda packages and performance of your system. You can allocate more threads if you prefer.
