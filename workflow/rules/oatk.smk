# This rule runs oatk to extract organelles reads from the hifi reads.
# get_oatk_outputs is used to dynamically decide the outputs of oatk

rule oatk:
    input:
        "results/reads/hifi/hifi.fastq.gz",
    output:
        gfa_mito = "results/oatk/oatk.asm.mito.gfa",
        fasta_mito = "results/oatk/oatk.asm.mito.ctg.fasta"
    threads: config["high"]["t"]
    log:
        "logs/oatk.log",
    conda:
        "../envs/oatk.yaml"
    resources:
        mem_mb=config["high"]["mem_mb"],  # access memory from config
    params:
        k=config["oatk"]["k"],
        c=config["oatk"]["c"],
        m=config["oatk"]["m"],
    shell:
        """
        /usr/bin/time -v oatk -k {params.k} -c {params.c} -t {threads} -m {params.m} {input} >> {log} 2>&1
        mv oatk.asm* results/oatk >> {log} 2>&1
        """


rule oatk_pltd:
    input:
        "results/reads/hifi/hifi.fastq.gz",
    output:
        gfa_mito = "results/oatk_pltd/oatk.asm.mito.gfa",
        fasta_mito = "results/oatk_pltd/oatk.asm.mito.ctg.fasta",
        gfa_pltd = "results/oatk_pltd/oatk.asm.pltd.gfa",
        fasta_pltd = "results/oatk_pltd/oatk.asm.pltd.ctg.fasta"
    threads: config["high"]["t"]
    log:
        "logs/oatk.log",
    conda:
        "../envs/oatk.yaml"
    resources:
        mem_mb=config["high"]["mem_mb"],  # access memory from config
    params:
        k=config["oatk"]["k"],
        c=config["oatk"]["c"],
        m=config["oatk"]["m"],
        optional_params=" ".join(
            f"{k} {v}" for k, v in config["oatk"]["optional_params"].items() if v
        ),
    shell:
        """
        /usr/bin/time -v oatk -k {params.k} -c {params.c} -t {threads} -m {params.m} {params.optional_params} {input} >> {log} 2>&1
        mv oatk.asm* results/oatk_pltd >> {log} 2>&1
        """