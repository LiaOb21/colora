# Define function for the inputs to rule_all
# this rule just defines the inputs to the final rule

def get_all_inputs():
    inputs = [
        "results/reads/hifi/hifi.fastq.gz",  # output of hifi_prep rule
        "results/nanoplot/NanoPlot-report.html",  # nanoplot report
        "results/kmc/out.hist",  # output of kmc rule
        "results/genomescope",  # output directory of genomescope
        "results/oatk/oatk.asm.mito.ctg.fasta",  # mitochondrion fasta from oatk
        "results/hifiasm/hifiasm.asm.p_ctg.gfa",  # draft hifiasm assembly.gfa
        "results/hifiasm/hifiasm.asm.p_ctg.fa",  # draft hifiasm assembly.fa 
        expand("results/fastp/{sample}_report_fastp.HiC.html", sample=samples), # fastp report for HiC reads pre-processing
        expand("results/arima_mapping_pipeline/RAW_DIR/{sample}_1.bam", sample=samples), # bwa mem bam file1
        expand("results/arima_mapping_pipeline/RAW_DIR/{sample}_2.bam", sample=samples), # bwa mem bam file2
        expand("results/arima_mapping_pipeline/FILT_DIR/{sample}_1.bam", sample=samples), # filtered bam file1
        expand("results/arima_mapping_pipeline/FILT_DIR/{sample}_2.bam", sample=samples), # filtered bam file2
        expand("results/arima_mapping_pipeline/TMP_DIR/{sample}.bam", sample=samples), # combined bam file
        expand("results/arima_mapping_pipeline/PAIR_DIR/{sample}.bam", sample=samples), # paired bam file
        expand("results/arima_mapping_pipeline/REP_DIR/{sample}_rep1.bam", sample=samples), # marked duplicates bam file
        expand("results/arima_mapping_pipeline/REP_DIR/metrics.{sample}_rep1.txt", sample=samples), # marked duplicates metrics file
        expand("results/arima_mapping_pipeline/REP_DIR/{sample}_rep1.bam.stats", sample=samples), # marked duplicates stats file
        expand("results/yahs/hifiasm_yahs_{sample}_scaffolds_final.fa", sample=samples) # final scaffolded assembly
    ]
    if config["oatk"]["optional_params"]["-p"]:
        inputs.append("results/oatk/oatk.asm.pltd.ctg.fasta")  # chloroplast fasta from oatk

    if config["include_fcsgx"] == True and config["include_purge_dups"] == False:
        inputs.append("results/ncbi_fcsgx/hifiasm.asm.p_ctg_clean.fa")  # decontaminated assembly
        inputs.extend([f"results/ncbi_fcsgx/hifiasm.asm.p_ctg_clean.fa.{ext}" for ext in ["sa", "amb", "ann", "bwt", "pac"]])  # bwa index outputs for decontaminated assembly
    elif config["include_purge_dups"] == True:
        inputs.append("results/purge_dups/hifiasm_p_purged.fa")  # purged assembly with purge_dups
        inputs.append("results/purge_dups_alt/hifiasm_a_purged.fa")  # purged alternate assembly with purge_dups2
        inputs.extend([f"results/purge_dups/hifiasm_p_purged.fa.{ext}" for ext in ["sa", "amb", "ann", "bwt", "pac"]])  # bwa index outputs of the purged primary assembly
    else:
        inputs.append("results/hifiasm/hifiasm.asm.p_ctg.fa")  # draft hifiasm assembly.fa 
        inputs.extend([f"results/hifiasm/hifiasm.asm.p_ctg.fa.{ext}" for ext in ["sa", "amb", "ann", "bwt", "pac"]])  # bwa index outputs for the draft hifiasm assembly

    return inputs


# create a function to read the inputs of purge_dups rule dynamically 
def get_purge_dups_inputs():
    if config["include_fcsgx"] == True:
        return "results/ncbi_fcsgx/hifiasm.asm.p_ctg_clean.fa"
    else:
        return "results/hifiasm/hifiasm.asm.p_ctg.fa"

# create a function to read the inputs of bwa_index rule dynamically 
# This will affect the "REF" input in the following rules:
# bwa_index
# bwa_mem
# two_read_bam_combiner
# yahs
def get_bwa_index_inputs():
    if config["include_purge_dups"] == True:
        return "results/purge_dups/hifiasm_p_purged.fa"
    elif config["include_fcsgx"] == True and config["include_purge_dups"] == False:
        return "results/ncbi_fcsgx/hifiasm.asm.p_ctg_clean.fa"    
    else:
        return "results/hifiasm/hifiasm.asm.p_ctg.fa"


# create a function to read the outputs of bwa_index rule dynamically 
# this rule is used also by bwa_mem to locate the idex files
def get_bwa_index_outputs():
    if config["include_purge_dups"] == True:
        return ["results/purge_dups/hifiasm_p_purged.fa.amb",
                "results/purge_dups/hifiasm_p_purged.fa.ann",
                "results/purge_dups/hifiasm_p_purged.fa.bwt",
                "results/purge_dups/hifiasm_p_purged.fa.pac",
                "results/purge_dups/hifiasm_p_purged.fa.sa"]
    elif config["include_fcsgx"] == True and config["include_purge_dups"] == False:
        return ["results/ncbi_fcsgx/hifiasm.asm.p_ctg_clean.fa.amb",
                "results/ncbi_fcsgx/hifiasm.asm.p_ctg_clean.fa.ann",
                "results/ncbi_fcsgx/hifiasm.asm.p_ctg_clean.fa.bwt",
                "results/ncbi_fcsgx/hifiasm.asm.p_ctg_clean.fa.pac",
                "results/ncbi_fcsgx/hifiasm.asm.p_ctg_clean.fa.sa"]
    else:
        return ["results/hifiasm/hifiasm.asm.p_ctg.fa.amb",
                "results/hifiasm/hifiasm.asm.p_ctg.fa.ann",
                "results/hifiasm/hifiasm.asm.p_ctg.fa.bwt",
                "results/hifiasm/hifiasm.asm.p_ctg.fa.pac",
                "results/hifiasm/hifiasm.asm.p_ctg.fa.sa"]
