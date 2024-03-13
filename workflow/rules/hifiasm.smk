# This rule uses hifiasm to assemble the genome from the hifi reads
# Following the instructions for obtaining primary/alternate assemblies: https://hifiasm.readthedocs.io/en/latest/pa-assembly.html#produce-primary-alternate-assemblies


rule hifiasm:
    input:
        "results/reads/hifi/hifi.fastq.gz",
    output:
        gfa="results/hifiasm/asm.primary.gfa",
        gfa_alt="results/hifiasm/asm.alternate.gfa",
        fasta="results/hifiasm/asm.primary.fa",
        fasta_alt="results/hifiasm/asm.alternate.fa",
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
        hifiasm {input} -t {threads} -o results/hifiasm/asm --primary {params.optional_params} >> {log} 2>&1
        mv asm.p_ctg.gfa {output.gfa}
        mv asm.a_ctg.gfa {output.gfa_alt}       
        awk -f scripts/gfa_to_fasta.awk < {output.gfa} > {output.fasta}
        awk -f scripts/gfa_to_fasta.awk < {output.gfa_alt} > {output.fasta_alt}
        """


rule hifiasm_het:
    input:
        "results/reads/hifi/hifi.fastq.gz",
    output:
        gfa="results/hifiasm/asm.hap1.gfa",
        gfa_alt="results/hifiasm/asm.hap2.gfa",
        fasta="results/hifiasm/asm.hap1.fa",
        fasta_alt="results/hifiasm/asm.hap2.fa",
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
        hifiasm {input} -t {threads} -o results/hifiasm/asm {params.optional_params} >> {log} 2>&1
        mv asm.hic.hap1.p_ctg.gfa {output.gfa}
        mv asm.hic.hap2.p_ctg.gfa {output.gfa_alt}
        awk -f scripts/gfa_to_fasta.awk < {output.gfa} > {output.fasta}
        awk -f scripts/gfa_to_fasta.awk < {output.gfa_alt} > {output.fasta_alt}
        """
