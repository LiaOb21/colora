# This rule includes the "Step 3.B: Add read group" and "Step 4: Mark duplicates" from the Arima Genomics mapping pipeline
# This pipeline is currently built for processing a single sample with one read1 and read2 FASTQ file.
# For more info: https://github.com/ArimaGenomics/mapping_pipeline
# colora doesn't contain the codes to handle technical and biological replicates of hic reads, refer to the original pipeline for that

rule picard:
    input:
        tmp_bam="results/arima_mapping_pipeline/TMP_DIR/{sample}.bam",
    output:
        bam_paired="results/arima_mapping_pipeline/PAIR_DIR/{sample}.bam",
        mark_dup="results/arima_mapping_pipeline/REP_DIR/{sample}_rep1.bam",
        metrics="results/arima_mapping_pipeline/REP_DIR/metrics.{sample}_rep1.txt",
        stats="results/arima_mapping_pipeline/REP_DIR/{sample}_rep1.bam.stats"
    log:
        "logs/{sample}_picard.log"
    conda:
        "../envs/arima_mapping_pipeline.yaml"
    shell:
        """
        PICARD_PATH=${{CONDA_PREFIX}}/share/picard-*/picard.jar
        mkdir -p results/arima_mapping_pipeline/PAIR_DIR/
        mkdir -p results/arima_mapping_pipeline/REP_DIR/
        java -Xmx4G -Djava.io.tmpdir=temp/ -jar ${PICARD} AddOrReplaceReadGroups INPUT={input.tmp_bam} OUTPUT={output.bam_paired} ID={wildcards.sample} LB={wildcards.sample} SM={wildcards.sample} PL=ILLUMINA PU=none 
        java -Xmx30G -XX:-UseGCOverheadLimit -Djava.io.tmpdir=temp/ -jar ${PICARD} MarkDuplicates INPUT={output.bam_paired} OUTPUT={output.mark_dup} METRICS_FILE={output.metrics} TMP_DIR=results/arima_mapping_pipeline/TMP_DIR ASSUME_SORTED=TRUE VALIDATION_STRINGENCY=LENIENT REMOVE_DUPLICATES=TRUE
        samtools index {output.mark_dup}
        perl scripts/get_stats.pl {output.mark_dup} > {output.stats}
        """