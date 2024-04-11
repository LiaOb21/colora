rule gfastats:
    input:
        "results/oatk/oatk.asm.mito.ctg.fasta"
    output:
        "results/gfastats/mito.stats"
    threads: config["low"]["t"]
    log:
        "logs/gfastats.log",
    resources:
        mem_mb=config["low"]["mem_mb"],  # access memory from config
    conda:
        "../envs/gfastats.yaml"
    shell:
        """
        (gfastats {input} > {output}) 2>> {log}
        """


rule gfastats_pltd:
    input:
        mito_in = "results/oatk_pltd/oatk.asm.mito.ctg.fasta",
        pltd_in = "results/oatk_pltd/oatk.asm.pltd.ctg.fasta"
    output:
        mito_out = "results/gfastats_pltd/mito.stats",
        pltd_out = "results/gfastats_pltd/pltd.stats"
    threads: config["low"]["t"]
    log:
        "logs/gfastats.log",
    resources:
        mem_mb=config["low"]["mem_mb"],  # access memory from config
    conda:
        "../envs/gfastats.yaml"
    shell:
        """
        (gfastats {input.mito_in} > {output.mito_out}) >> {log} 2>&1
        (gfastats {input.pltd_in} > {output.pltd_out}) >> {log} 2>&1
        """