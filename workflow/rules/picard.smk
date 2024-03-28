# This rule includes the "Step 3.B: Add read group" and "Step 4: Mark duplicates" from the Arima Genomics mapping pipeline
# This pipeline is currently built for processing a single sample with one read1 and read2 FASTQ file.
# For more info: https://github.com/ArimaGenomics/mapping_pipeline
# colora doesn't contain the codes to handle technical and biological replicates of hic reads, refer to the original pipeline for that
# The final files will be in REP_DIR to be consistent with the original pipeline but they will be marked with "final" and not with "rep1"


rule picard:
    input:
        tmp_bam="results/arima_mapping_pipeline_{hap}/TMP_DIR/hic_vs_contigs_filt_paired.bam",
    output:
        bam_add_read_group="results/arima_mapping_pipeline_{hap}/PAIR_DIR/paired_add_read_group.bam",
        mark_dup="results/arima_mapping_pipeline_{hap}/REP_DIR/paired_mark_dups_final.bam",
        metrics="results/arima_mapping_pipeline_{hap}/REP_DIR/metrics.final.txt",
        stats="results/arima_mapping_pipeline_{hap}/REP_DIR/final.bam.stats",
    params:
        hic_sample_name=sample
    log:
        "logs/picard_{hap}.log",
    resources:
        mem_mb=config['arima']['mem_mb'],  # access memory from config
    conda:
        "../envs/arima_mapping_pipeline.yaml"
    shell:
        """
        PICARD=${{CONDA_PREFIX}}/share/picard-*/picard.jar
        java -Xmx4G -Djava.io.tmpdir=temp/ -jar ${{PICARD}} AddOrReplaceReadGroups INPUT={input.tmp_bam} OUTPUT={output.bam_add_read_group} ID={params.hic_sample_name} LB={params.hic_sample_name} SM={params.hic_sample_name} PL=ILLUMINA PU=none 
        java -Xmx30G -XX:-UseGCOverheadLimit -Djava.io.tmpdir=temp/ -jar ${{PICARD}} MarkDuplicates INPUT={output.bam_add_read_group} OUTPUT={output.mark_dup} METRICS_FILE={output.metrics} TMP_DIR=results/arima_mapping_pipeline_{wildcards.hap}/TMP_DIR ASSUME_SORTED=TRUE VALIDATION_STRINGENCY=LENIENT REMOVE_DUPLICATES=TRUE
        samtools index {output.mark_dup}
        perl scripts/get_stats.pl {output.mark_dup} > {output.stats}
        """
