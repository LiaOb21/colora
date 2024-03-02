# This rule uses hifiasm to assemble the genome from the hifi reads
# Following the instructions for obtaining primary/alternate assemblies: https://hifiasm.readthedocs.io/en/latest/pa-assembly.html#produce-primary-alternate-assemblies


rule hifiasm:
    input:
        "results/reads/hifi/hifi.fastq.gz",
    output:
        gfa="results/hifiasm/hifiasm.asm.p_ctg.gfa",
        gfa_alt="results/hifiasm/hifiasm.asm.a_ctg.gfa",
        fasta="results/hifiasm/hifiasm.asm.p_ctg.fa",
        fasta_alt="results/hifiasm/hifiasm.asm.a_ctg.fa",
    threads: config["hifiasm"]["t"]
    log:
        "logs/hifiasm.log",
    resources:
        mem_mb=config['hifiasm']['mem_mb'],  # access memory from config
    conda:
        "../envs/hifiasm.yaml"
    params:
        optional_params=" ".join(
            f"{k} {v}" for k, v in config["hifiasm"]["optional_params"].items() if v
        ),
    shell:
        """
        hifiasm {input} -t {threads} -o results/hifiasm/hifiasm.asm --primary {params.optional_params} >> {log} 2>&1
        awk -f scripts/gfa_to_fasta.awk < {output.gfa} > {output.fasta}
        awk -f scripts/gfa_to_fasta.awk < {output.gfa_alt} > {output.fasta_alt}
        """
