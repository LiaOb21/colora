# This rule includes the "Step 3A: Pair reads & mapping quality filter" from the Arima Genomics mapping pipeline
# This pipeline is currently built for processing a single sample with one read1 and read2 FASTQ file.
# For more info: https://github.com/ArimaGenomics/mapping_pipeline
# colora doesn't contain the codes to handle technical and biological replicates of hic reads, refer to the original pipeline for that



rule two_read_bam_combiner:
    input:
        bam1_filt="results/arima_mapping_pipeline_{hap}/FILT_DIR/hic_vs_contigs_filt_1.bam",
        bam2_filt="results/arima_mapping_pipeline_{hap}/FILT_DIR/hic_vs_contigs_filt_2.bam",
        REF="results/bwa_index_{hap}/asm.fa",
    output:
        tmp_bam="results/arima_mapping_pipeline_{hap}/TMP_DIR/hic_vs_contigs_filt_paired.bam",
    threads: config["medium"]["t"]
    params:
        MAPQ_FILTER=config["arima"]["MAPQ_FILTER"],
    log:
        "logs/two_read_bam_combiner_{hap}.log",
    resources:
        mem_mb=config["medium"]["mem_mb"],  # access memory from config
    conda:
        "../envs/arima_mapping_pipeline.yaml"
    shell:
        """
        (samtools faidx {input.REF}) >> {log} 2>&1
        (perl scripts/two_read_bam_combiner.pl {input.bam1_filt} {input.bam2_filt} samtools {params.MAPQ_FILTER} | samtools view -bS -t {input.REF}.fai - | samtools sort -@ {threads} -o {output.tmp_bam} - ) >> {log} 2>&1
        """
