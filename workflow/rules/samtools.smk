# This rule uses samtools and shell commands to extract the name of the fastq reads aligning to organelles fasta for their subsequent removal

rule run_samtools:
    input:
        mito = "results/minimap2/aln.mito.sam",
        pltd = "results/minimap2/aln.pltd.sam"
    output:
        mito_bam = "results/samtools/mito.sorted.bam",
        pltd_bam = "results/samtools/pltd.sorted.bam",
        mito_reads = "results/samtools/mapped.mito.reads.txt",
        pltd_reads = "results/samtools/mapped.pltd.reads.txt",
        organelles = "results/samtools/organelles.reads.txt"
    conda:
        "../envs/samtools.yaml"
    shell:
        """
        samtools view -Sb {input.mito} | samtools sort -o {output.mito_bam}
        samtools view -Sb {input.pltd} | samtools sort -o {output.pltd_bam}
        samtools view -F 4 {output.mito_bam} | cut -f1 > {output.mito_reads}
        samtools view -F 4 {output.pltd_bam} | cut -f1 > {output.pltd_reads}
        cat {output.mito_reads} {output.pltd_reads} | sort | uniq > {output.organelles}
        """