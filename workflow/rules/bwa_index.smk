# This rule includes the "Step 0: Index reference" from the Arima Genomics mapping pipeline
# This pipeline is currently built for processing a single sample with one read1 and read2 FASTQ file.
# For more info: https://github.com/ArimaGenomics/mapping_pipeline
# colora doesn't contain the codes to handle technical and biological replicates of hic reads, refer to the original pipeline for that


# include common.smk to use get_bwa_index_inputs and get_bwa_index_outputs functions
include: "common.smk"


rule bwa_index:
    input:
        REF=get_bwa_index_inputs(),
    output:
        get_bwa_index_outputs(),
    log:
        "logs/bwa_index.log",
    conda:
        "../envs/arima_mapping_pipeline.yaml"
    shell:
        """
        bwa index -a bwtsw {input.REF} >> {log} 2>&1
        """
