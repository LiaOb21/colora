# This rule uses seqkit to filter organelle reads from the initial fastq files

rule run_seqkit:
    input:
        organelles = "results/samtools/organelles.reads.txt",
        fastq = "results/reads/hifi/hifi.fastq.gz"
    output:
        organelles_filtered = "results/reads/hifi/filtered_orgs_reads.fastq.gz"
    conda:
        "../envs/seqkit.yaml"
    shell:
        """
        seqkit grep -f {input.organelles} -v {input.fastq} | gzip > {output.organelles_filtered}
        """