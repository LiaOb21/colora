# This rule includes the "Step 0: Index reference" from the Arima Genomics mapping pipeline
# This pipeline is currently built for processing a single sample with one read1 and read2 FASTQ file.
# For more info: https://github.com/ArimaGenomics/mapping_pipeline
# colora doesn't contain the codes to handle technical and biological replicates of hic reads, refer to the original pipeline for that

rule bwa_index:
    input:
        get_bwa_index_inputs
    output:
        dir = directory("results/bwa_index_{hap}"),
        asm = "results/bwa_index_{hap}/asm.fa"
    log:
        "logs/bwa_index_{hap}.log",
    resources:
        mem_mb=config["low"]["mem_mb"],  # access memory from config
    conda:
        "../envs/arima_mapping_pipeline.yaml"
    shell:
        """
        ln -srn {input} {output.asm}
        /usr/bin/time -v bwa index -a bwtsw {output.asm} >> {log} 2>&1
        """