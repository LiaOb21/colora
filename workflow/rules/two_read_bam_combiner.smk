# This rule includes the "Step 3A: Pair reads & mapping quality filter" from the Arima Genomics mapping pipeline
# This pipeline is currently built for processing a single sample with one read1 and read2 FASTQ file.
# For more info: https://github.com/ArimaGenomics/mapping_pipeline
# colora doesn't contain the codes to handle technical and biological replicates of hic reads, refer to the original pipeline for that

rule two_read_bam_combiner:
    input:
        bam1_filt="results/arima_mapping_pipeline/FILT_DIR/{sample}_1.bam",
        bam2_filt="results/arima_mapping_pipeline/FILT_DIR/{sample}_2.bam",
        REF = "results/purge_dups/hifiasm_p_purged.fa"
    output:
        tmp_bam="results/arima_mapping_pipeline/TMP_DIR/{sample}.bam"
    log:
        "logs/{sample}_two_read_bam_combiner.log"
    conda:
        "../envs/arima_mapping_pipeline.yaml"
    shell:
        """
        perl scripts/two_read_bam_combiner.pl {input.bam1_filt} {input.bam2_filt} samtools {config[arima][MAPQ_FILTER]} | samtools view -bS -t {input.REF}.fai - | samtools sort -@ {config[arima][CPU]} -o {output.tmp_bam} - >> {log} 2>&1
        """