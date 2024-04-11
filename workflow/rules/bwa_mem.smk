# This rule includes the "Step 1.A: FASTQ to BAM (1st)" and "Step 1.B: FASTQ to BAM (2nd)" from the Arima Genomics mapping pipeline
# This pipeline is currently built for processing a single sample with one read1 and read2 FASTQ file.
# For more info: https://github.com/ArimaGenomics/mapping_pipeline
# colora doesn't contain the codes to handle technical and biological replicates of hic reads, refer to the original pipeline for that



rule bwa_mem:
    input:
        REF="results/bwa_index_{hap}/asm.fa",
        forward_hic="results/fastp/hic_trim_1.fastq.gz",
        reverse_hic="results/fastp/hic_trim_2.fastq.gz",
    output:
        bam1="results/arima_mapping_pipeline_{hap}/RAW_DIR/hic_vs_contigs_1.bam",
        bam2="results/arima_mapping_pipeline_{hap}/RAW_DIR/hic_vs_contigs_2.bam",
    threads: config["arima"]["CPU"]
    log:
        "logs/bwa_mem_{hap}.log",
    resources:
        mem_mb=config['arima']['mem_mb'],  # access memory from config    
    conda:
        "../envs/arima_mapping_pipeline.yaml"
    shell:
        """
        /usr/bin/time -v sh -c 'bwa mem -M -t {threads} {input.REF} {input.forward_hic} | samtools view -@ {threads} -Sb - > {output.bam1}' >> {log} 2>&1
        /usr/bin/time -v sh -c 'bwa mem -M -t {threads} {input.REF} {input.reverse_hic} | samtools view -@ {threads} -Sb - > {output.bam2}' >> {log} 2>&1
        """
