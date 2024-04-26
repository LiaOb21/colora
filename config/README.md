
***Please, read carefully.*** :snake:

# Setting for memory and threads

The rules have been divided in high demanding, medium demanding and low demanding based on testing with different species.

In this part of the `config.yaml` we can set these values.

### 1. Example of suitable values for barley in our testing (genome size ~4.2 Gbp, 75% of repeats):

```
# Set memory and threads for high demanding rules
high:
  mem_mb: 409600 # memory in MB
  t: 50 # number of threads

# Set memory and threads for medium demanding rules
medium:
  mem_mb: 204800 # memory in MB
  t: 20 # number of threads

# Set memory and threads for low demanding rules
low:
  mem_mb: 51200 # memory in MB
  t: 8 # number of threads
```

# Paths to reads

Below is an example of how to set the paths to hifi and hic reads. It's showing the case in which hifi and hic have been placed under `resources/raw_hifi` and `resources/raw_hic` directories into `colora` directory itself, respectively. It's also possible to specify different paths if you prefer, but please notice that if you are not using subdirectories of `colora` itself you should use absolute paths.

```
# Path to hifi reads
hifi_path: "resources/raw_hifi/"

# Path to hic reads
hic_path: "resources/raw_hic/"
```

# KMC and genomescope2

KMC produces the k-mer database that is then used by genomescope2 to have an estimate of the genome sie and other statistics. Therefore, parameters set for KMC are going to influence the behaviour of genomescope2 as well.

The k-mer size, if unsure, can be estimated using a script provided by Merqury, called `best_k.sh`, available [here](https://github.com/marbl/merqury/blob/master/best_k.sh).

For the parameter `cs` you must pay attention to your genome size and repetitiveness, as this can lead to wrong estimation of genome size. See [genomescope2 FAQs](https://github.com/tbenavi1/genomescope2.0?tab=readme-ov-file#frequently-asked-questions-faq).

### 1. Example of suitable parameters for barley in our testing (genome size ~4.2 Gbp, 75% of repeats):

```
# Customisable parameters for kmc
kmc:
  k: 21 # kmer size, it will be the same used for genomescope2
  ci: 1 # exclude k-mers occurring less than <value> times (default: 2)
  cs: 1000000 #maximal value of a counter (default: 255)

# Customisable parameters for kmc_tools transform
kmc_tools:
  cx: 1000000 # exclude k-mers occurring more of than <value> times
```

# oatk

See usage section for instructions on how to obtain the oatk database of interest.

Note that the `c` parameter, which specifies coverage threshold, can be set as 5-10 times the value of nuclear sequence coverage.

### 1. Example of suitable parameters for barley in our testing (hifi coverage ~30x):

```
# Customisable parameters for oatk
oatk:
  k: 1001 # kmer size [1001]
  c: 150 #  minimum kmer coverage [3]
  m: "resources/oatkDB/embryophyta_mito.fam" # mitochondria gene annotation HMM profile database [NULL]
  optional_params: 
    "-p": "resources/oatkDB/embryophyta_pltd.fam" # to use for species that have a plastid db
```

# fastp

Fastp is used for the preprocessing of the hic reads. Adapters are always removed, while it is suggested to use the optional parameters in the `config.yaml` when the data are obtained with Arima Hi-C library prep kit.

### 1. Example of parameters to be used with data obtained with Arima Hi-C library prep kit:

```
# Customisable parameters for fastp
fastp:
  optional_params:
    "--cut_front": True # to use only with Arima Hi-C library prep kit generated data
    "--cut_front_window_size": "5" # to use only with Arima Hi-C library prep kit generated data
```

### 2. Example of parameters to be used with other library prep kit:

```
# Customisable parameters for fastp
fastp:
  optional_params:
    "--cut_front": False # to use only with Arima Hi-C library prep kit generated data
    "--cut_front_window_size": "" # to use only with Arima Hi-C library prep kit generated data
```

# hifiasm

Hifiasm is the assembler used by colora. It is possible to obtain a primary assembly or a phased assembly with colora. You can choose this by setting `phased_assembly: False` or `phased_assembly: True` respectively in the `config.yaml`. In the case of phased assembly, Hi-C reads are integrated in the assembly process, after being pre-processed with fastp to remove adapters (and 5' ends in the case of Arima library prep kit).

All the other parameters are optional, and you can refer to [hifiasm documentation](https://hifiasm.readthedocs.io/en/latest/index.html) to understand what is optimal for your case.

If you also have ONT reads for your sample, you may consider using them in  the assembly process. Remember that you should previously trim them (if needed) and join them in a single file (if they are in multiple files).

### 1. Example of primary assembly (hifiasm uses only HiFi reads in this case):

```
# Customisable parameters for hifiasm
hifiasm:
  phased_assembly: False # set to true if you want to obtain a phased assembly
  optional_params:
    "-f": "" # used for small datasets
    "-l": "" # purge level. 0: no purging; 1: light; 2/3: aggressive [0 for trio; 3 for unzip]
    "--ul": "" # use this if you have also ont data yu want to integrate in your assembly
```

### 2. Example of primary assembly using HiFi and ONT reads:

```
# Customisable parameters for hifiasm
hifiasm:
  phased_assembly: False # set to true if you want to obtain a phased assembly
  optional_params:
    "-f": "" # used for small datasets
    "-l": "" # purge level. 0: no purging; 1: light; 2/3: aggressive [0 for trio; 3 for unzip]
    "--ul": "resources/raw_ont/ont_reads.fastq.gz" # use this if you have also ont data you want to integrate in your assembly
```

### 3. Example of phased assembly with extra optional parameters (in this case HiFi and Hi-C reads are used):

```
# Customisable parameters for hifiasm
hifiasm:
  phased_assembly: True # set to true if you want to obtain a phased assembly
  optional_params:
    "--hom-cov": "181"
    "-f": "" # used for small datasets
    "-l": "" # purge level. 0: no purging; 1: light; 2/3: aggressive [0 for trio; 3 for unzip]
    "--ul": "" # use this if you have also ont data yu want to integrate in your assembly
```

### 4. Example of phased assembly with ONT integration (in this case HiFi, Hi-C and ONT reads are used):

```
# Customisable parameters for hifiasm
hifiasm:
  phased_assembly: True # set to true if you want to obtain a phased assembly
  optional_params:
    "-f": "" # used for small datasets
    "-l": "" # purge level. 0: no purging; 1: light; 2/3: aggressive [0 for trio; 3 for unzip]
    "--ul": "resources/raw_ont/ont_reads.fastq.gz" # use this if you have also ont data yu want to integrate in your assembly
```

# FCS-gx

FCS-gx is the tool that performs the decontamination of the assembly. However, this rule can be skipped setting `include_fcsgx: False`. This choice is made given the fact that the database that FCS-gx uses is ~500GB, and it needs a big RAM, and we want colora to be as portable as possible. Skipping this rule it's possible to run colora even in a small system.

### 1. Example of suitable parameters for barley:

```
#Set this to False if you want to skip the fcsgx step:
include_fcsgx: True #inlcude this rule only if you have preiously downloaded the database (recommended to run fcsgx only on a HPC. It requires around 500 GB of space on your disk and a large RAM)

# Customisable parameters for fcsgx
fcsgx:
  ncbi_tax_id: 4513
  path_to_gx_db: "resources/gxdb"
```

### 2. Example of suitable parameters for Rhizophagus irregularis:

```
#Set this to False if you want to skip the fcsgx step:
include_fcsgx: True #inlcude this rule only if you have preiously downloaded the database (recommended to run fcsgx only on a HPC. It requires around 500 GB of space on your disk and a large RAM)

# Customisable parameters for fcsgx
fcsgx:
  ncbi_tax_id: 588596
  path_to_gx_db: "resources/gxdb"
```

# purge_dups

Purge_dups is a tool used in most of the assembly pipelines to purge haplotigs and overlaps in an assembly based on read depth. In `colora`, the rule used to run this tool is run only when `phased_assembly: False` is set for hifiasm, it is optional and can be skipped setting `include_purge_dups: False`. This is because in some of our projects we found that purge_dups was excluding biological relevant information from the assembly, giving worse busco scores compared to the assembly obtained skipping this rule. However, it may be beneficial to try several settings (including purge_dups or not) to see what combinations is the best for your case. 

When `phased_assembly: True` is set for hifiasm, it is mandatory to set `include_purge_dups: False`. This is because the purging step is not necessary when we obtain two assembly for two different haplotypes. 

When `phased_assembly: False` is set for hifiasm, you can choose to set `include_purge_dups: True` or `include_purge_dups: False`. Again, we recommend to try both ways to see what strategy gives better results.

# Arima genomics mapping pipeline

This pipeline has been integrated in colora, translating it from the [original bash pipeline](https://github.com/ArimaGenomics/mapping_pipeline) to a snakemake version. 

The only changes we made compared to the original pipeline are the following:

- Remove the PREFIX line and the option -p $PREFIX from the bwa command (step 0 of the original pipeline)
- add -M flag in bwa mem command (step 1.A and 1.B of the original pipeline)
- pipeline split in several rules

Moreover, colora doesn't handle technical and biological replicates of Hi-C reads, so you must refer to the original pipeline for that. 

The preprocessing of the Hi-C reads is made with fastp (see above).

Samples can be processed only one by one with this pipeline.

If you set `phased_assembly: True` for hifiasm, the two assemblies for the different haplotypes are processed in parallel.

The only parameter you can custumise for the arima genomics mapping pipeline is the mapping quality filter (10 is the default):

### 1. Example of default mapping quality filter for arima mapping pipeline

```
# Customisable parameters for arima mapping pipeline:
arima:
  MAPQ_FILTER: 10
```

# yahs

YaHS is the tool that we use for performing the scaffolding on the contig level assembly using Hi-C reads. 

The only customisable parameter for yahs in `colora` is the restriction enzyme(s), but it is an optional parameter and you can choose to set it or not:

```
# Customisable parameters for yahs
yahs:
  optional_params: 
    "-e": "" # you can specify the restriction enzyme(s) used by the Hi-C experiment
```

N.B. From [YaHS documentation](https://github.com/c-zhou/yahs): 

> With -e option, you can specify the restriction enzyme(s) used by the Hi-C experiment. For example, GATC for the DpnII restriction enzyme used by the Dovetail Hi-C Kit; GATC,GANTC and GATC,GANTC,CTNAG,TTAA for Arima genomics 2-enzyme and 4-enzyme protocol, respectively. Sometimes, the specification of enzymes may not change the scaffolding result very much if not make it worse, especially when the base quality of the assembly is not very good, e.g., assembies constructed from noisy long reads.


# Quast

Quast is used to perform the quality check of the assembly along all the workflow. It processes:

- contig level assembly(ies) produced by hifiasm
- decontaminated contig level assembly(ies) produced by fcs-gx (if included)
- purged contig level primary assembly produced by purge_dups (if included)
- scaffolded assembly(ies) produced by yahs

You can optionally give as input the reference genome for the species under study to compare you assembly with.

All the parameters that you can customise for quast are optional, and you can refer to [quast documentation](https://quast.sourceforge.net/docs/manual.html) to see which ones are suitable for your case.

### 1. Example
### 2. Example

# Busco
