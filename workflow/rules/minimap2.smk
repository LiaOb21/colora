# This rule uses minimap2 to align the fastq reads to the organelles fasta

rule run_minimap2:
    input:
        mito = "results/oatk/oatk.asm.mito.ctg.fasta",
        pltd = "results/oatk/oatk.asm.pltd.ctg.fasta",
        fastq = "results/reads/hifi/hifi.fastq.gz"
    output:
        mito_sam = "results/minimap2/aln.mito.sam",
        pltd_sam = "results/minimap2/aln.pltd.sam",
    conda:
        "../envs/minimap2.yaml"
    shell:
        """
        minimap2 -ax map-hifi {input.mito} -t {config[minimap2][t]} {input.fastq} > {output.mito_sam}
        minimap2 -ax map-hifi {input.pltd} -t {config[minimap2][t]} {input.fastq} > {output.pltd_sam}
        """