# This rule includes the "Step 1.A: FASTQ to BAM (1st)" and "Step 1.B: FASTQ to BAM (2nd)" from the Arima Genomics mapping pipeline
# This pipeline is currently built for processing a single sample with one read1 and read2 FASTQ file.
# For more info: https://github.com/ArimaGenomics/mapping_pipeline
# colora doesn't contain the codes to handle technical and biological replicates of hic reads, refer to the original pipeline for that

rule bwa_mem:
    input:
        REF = "results/purge_dups/hifiasm_p_purged.fa",
        index = rules.bwa_index.output[0],
        forward_hic = "results/fastp/{sample}_trim_1.fastq.gz",
        reverse_hic = "results/fastp/{sample}_trim_2.fastq.gz"
    output:
        bam1="results/arima_mapping_pipeline/RAW_DIR/{sample}_1.bam",
        bam2="results/arima_mapping_pipeline/RAW_DIR/{sample}_2.bam"
    log:
        "logs/{sample}_bwa_mem.log"
    conda:
        "../envs/arima_mapping_pipeline.yaml"
    shell:
        """
        mkdir -p results/arima_mapping_pipeline/RAW_DIR
        bwa mem -M -t {config[arima][CPU]} {input.REF} {input.forward_hic} | samtools view -@ {config[arima][CPU]} -Sb - > {output.bam1} 
        bwa mem -M -t {config[arima][CPU]} {input.REF} {input.reverse_hic} | samtools view -@ {config[arima][CPU]} -Sb - > {output.bam2} 
        """