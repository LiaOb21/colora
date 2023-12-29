# This rule uses hifiasm to assemble the genome after filtering out organelles reads
# Following the instructions for assembling heterozygous genomes from the hifiasm github

GFA_TO_FASTA = r"""/^S/{print ">"$2;print $3}"""

rule run_hifiasm:
    input:
        "results/reads/hifi/hifi.fastq.gz",
    output:
        gfa="results/hifiasm/hifiasm.asm.bp.p_ctg.gfa",
        fasta="results/hifiasm/hifiasm.asm.bp.p_ctg.fa"
    log:
        "logs/hifiasm.log"
    conda:
        "../envs/hifiasm.yaml"
    params:
        optional_params=" ".join(
            f"{k} {v}" for k, v in config["hifiasm"]["optional_params"].items() if v
        )
    shell:
        """
        hifiasm {input} -t {config[hifiasm][t]} -o results/hifiasm/hifiasm.asm {params.optional_params} 
        awk {GFA_TO_FASTA:q} {output.gfa} > {output.fasta}
        """
