# This rule uses hifiasm to assemble the genome from the hifi reads
# Following the instructions for obtaining primary/alternate assemblies: https://hifiasm.readthedocs.io/en/latest/pa-assembly.html#produce-primary-alternate-assemblies


rule hifiasm:
    input:
        reads = "results/reads/hifi/hifi.fastq.gz",
    output:
        gfa="results/hifiasm/asm.primary.gfa",
        gfa_alt="results/hifiasm/asm.alternate.gfa",
        fasta="results/hifiasm/asm.primary.fa",
        fasta_alt="results/hifiasm/asm.alternate.fa",
        link_primary = "results/assemblies/asm_primary.fa",
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
        hifiasm {input.reads} -t {threads} -o results/hifiasm/asm --primary {params.optional_params} >> {log} 2>&1
        mv results/hifiasm/asm.p_ctg.gfa {output.gfa}
        mv results/hifiasm/asm.a_ctg.gfa {output.gfa_alt}       
        awk -f scripts/gfa_to_fasta.awk < {output.gfa} > {output.fasta}
        awk -f scripts/gfa_to_fasta.awk < {output.gfa_alt} > {output.fasta_alt}

        # all the assemblies produced by the workflow will be symlinked to results/assemblies

        ln -srn {output.fasta} {output.link_primary}
        """


rule hifiasm_het:
    input:
        reads = "results/reads/hifi/hifi.fastq.gz",
        hic_forward = "results/fastp/hic_trim_1.fastq.gz",
        hic_reverse = "results/fastp/hic_trim_2.fastq.gz"
    output:
        gfa_hap1="results/hifiasm/asm.hap1.gfa",
        gfa_hap2="results/hifiasm/asm.hap2.gfa",
        fasta_hap1="results/hifiasm/asm.hap1.fa",
        fasta_hap2="results/hifiasm/asm.hap2.fa",
        link_hap1 = "results/assemblies/asm_hap1.fa",
        link_hap2 = "results/assemblies/asm_hap2.fa"
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
        hifiasm {input.reads} -t {threads} -o results/hifiasm/asm --h1 {input.hic_forward} --h2 {input.hic_reverse} {params.optional_params} >> {log} 2>&1
        mv results/hifiasm/asm.hic.hap1.p_ctg.gfa {output.gfa_hap1}
        mv results/hifiasm/asm.hic.hap2.p_ctg.gfa {output.gfa_hap2}
        awk -f scripts/gfa_to_fasta.awk < {output.gfa_hap1} > {output.fasta_hap1}
        awk -f scripts/gfa_to_fasta.awk < {output.gfa_hap2} > {output.fasta_hap2}

        # all the assemblies produced by the workflow will be symlinked to results/assemblies

        ln -srn {output.fasta_hap1} {output.link_hap1}
        ln -srn {output.fasta_hap2} {output.link_hap2}       
        """
