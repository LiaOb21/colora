# This rule includes the "Step 2.A: Filter 5' end (1st)" and "Step 2.B: Filter 5' end (2nd)" from the Arima Genomics mapping pipeline
# This pipeline is currently built for processing a single sample with one read1 and read2 FASTQ file.
# For more info: https://github.com/ArimaGenomics/mapping_pipeline
# colora doesn't contain the codes to handle technical and biological replicates of hic reads, refer to the original pipeline for that

rule fiter_five_end:
    input:
        bam1="results/arima_mapping_pipeline/RAW_DIR/{sample}_1.bam",
        bam2="results/arima_mapping_pipeline/RAW_DIR/{sample}_2.bam"
    output:
        bam1_filt="results/arima_mapping_pipeline/FILT_DIR/{sample}_1.bam",
        bam2_filt="results/arima_mapping_pipeline/FILT_DIR/{sample}_2.bam"
    log:
        "logs/{sample}_fiter_five_end.log"
    conda:
        "../envs/arima_mapping_pipeline.yaml"
    shell:
        """
        samtools view -h {input.bam1} | perl scripts/filter_five_end.pl | samtools view -Sb - > {output.bam1_filt} >> {log} 2>&1
        samtools view -h {input.bam2} | perl scripts/filter_five_end.pl | samtools view -Sb - > {output.bam2_filt} >> {log} 2>&1
        """