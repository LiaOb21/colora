# This rule includes the "Step 0: Index reference" from the Arima Genomics mapping pipeline
# This pipeline is currently built for processing a single sample with one read1 and read2 FASTQ file.
# For more info: https://github.com/ArimaGenomics/mapping_pipeline
# colora doesn't contain the codes to handle technical and biological replicates of hic reads, refer to the original pipeline for that


# create a function to read the inputs of bwa_index rule dynamically 
def get_bwa_index_inputs():
    if config["include_purge_dups"] == 'True':
        return "results/purge_dups/hifiasm_p_purged.fa"
    else:
        return "results/hifiasm/hifiasm.asm.p_ctg.fa"


# create a function to read the outputs of bwa_index rule dynamically 
def get_bwa_index_outputs():
    if config["include_purge_dups"] == 'True':
        return ["results/purge_dups/hifiasm_p_purged.fa.amb",
                "results/purge_dups/hifiasm_p_purged.fa.ann",
                "results/purge_dups/hifiasm_p_purged.fa.bwt",
                "results/purge_dups/hifiasm_p_purged.fa.pac",
                "results/purge_dups/hifiasm_p_purged.fa.sa"]
    else:
        return ["results/hifiasm/hifiasm.asm.p_ctg.fa.amb",
                "results/hifiasm/hifiasm.asm.p_ctg.fa.ann",
                "results/hifiasm/hifiasm.asm.p_ctg.fa.bwt",
                "results/hifiasm/hifiasm.asm.p_ctg.fa.pac",
                "results/hifiasm/hifiasm.asm.p_ctg.fa.sa"]

rule bwa_index:
    input:
        REF = get_bwa_index_inputs()
    output:
        get_bwa_index_outputs()
    log:
        "logs/bwa_index.log"
    conda:
        "../envs/arima_mapping_pipeline.yaml"
    shell:
        """
        bwa index -a bwtsw {input.REF} >> {log} 2>&1
        """