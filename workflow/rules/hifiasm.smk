# This rule uses hifiasm to assemble the genome after filtering out organelles reads
# Following the instructions for assembling heterozygous genomes from the hifiasm github


rule run_hifiasm:
    input:
        organelles_filtered="results/reads/hifi/filtered_orgs_reads.fastq.gz",
    output:
        gfa="results/assemblies/hifiasm/hifiasm.asm.bp.p_ctg.gfa",
        fasta="results/assemblies/hifiasm/hifiasm.asm.bp.p_ctg.fa",
    conda:
        "../envs/hifiasm.yaml"
    params:
        optional_params=" ".join(
            f"{k} {v}" for k, v in config["hifiasm"]["optional_params"].items() if v
        ),
    shell:
        """
        hifiasm {input.organelles_filtered} -t {config[hifiasm][t]} -o results/assemblies/hifiasm/hifiasm.asm {params.optional_params} 
        awk '/^S/{{print ">"$$2;print $$3}}' {output.gfa} > {output.fasta}
        """
