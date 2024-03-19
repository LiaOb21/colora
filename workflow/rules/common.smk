# Define function for the inputs to rule_all
# this rule just defines the inputs to the final rule

import os

wildcard_constraints:
    hap = "primary|hap1|hap2"

def get_all_inputs(wc=None):
    if diploid_mode: 
        hap= ["hap1", "hap2"]
    else:
        hap = ["primary"]
    if not samples:
        # Show the user a meaningful error message
        logger.error(f"No files matching {samples_in_pattern} were found in"
                     f" {config.get('hic_path')}. Please check your config file.")
        logger.error(f"Exiting.")
        sys.exit(1)
    if diploid_mode and config["include_purge_dups"] == True:
        raise ValueError("--h1 and --h2 are set in the hifiasm rule and this means that you want to obtain a phased assembly. In this case you must skip purge_dups, so please set `include_purge_dups: False` in your config.yaml.")

    inputs = [
        "results/nanoplot/NanoPlot-report.html",  # nanoplot report
        "results/kmc/out.hist",  # output of kmc rule
        "results/genomescope",  # output directory of genomescope
    ]

    # Get the OATK outputs
    oatk_outputs = get_oatk_outputs()

    # If the OATK outputs is a string, convert it to a list
    if isinstance(oatk_outputs, str):
        oatk_outputs = [oatk_outputs]

    # Extend the list of inputs with the OATK outputs
    inputs.extend(oatk_outputs)

    inputs.extend(expand(
            "results/yahs_{hap}/asm_yahs_{sample}_scaffolds_final.fa", sample=samples, hap=hap
        ))  # final scaffolded assembly

    # If not in diploid mode, add asm.alternate.fa as an input
    if not diploid_mode:
        inputs.append("results/purge_dups_alt/asm.alternate_purged.fa")

    return inputs


# create a function to read the outputs of oatk rule dynamically (if the plastid database is provided or not)
def get_oatk_outputs():
    if config["oatk"]["optional_params"]["-p"]:
        return [
            "results/oatk/oatk.asm.mito.ctg.fasta",
            "results/oatk/oatk.asm.pltd.ctg.fasta",
        ]
    else:
        return "results/oatk/oatk.asm.mito.ctg.fasta"

# create a function to read the inputs of purge_dups rule dynamically
def get_purge_dups_inputs():
    if config["include_fcsgx"] == True:
        return "results/ncbi_fcsgx_primary/asm_clean.fa"
    else:
        return "results/hifiasm/asm.primary.fa"


# create a function to read the inputs of bwa_index rule dynamically
# This will affect the "REF" input in the following rules:
# bwa_index
# bwa_mem
# two_read_bam_combiner
# yahs

def get_bwa_index_inputs(wildcards):
    hap = wildcards.hap
    if config["include_purge_dups"] == True:
        return "results/purge_dups/asm.{hap}_purged.fa"
    elif config["include_fcsgx"] == True:
        return f"results/ncbi_fcsgx_{hap}/asm_clean.fa"
    else:
        return f"results/hifiasm/asm.{hap}.fa"


