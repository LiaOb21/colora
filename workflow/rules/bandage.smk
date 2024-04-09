rule bandage:
    input:
        "results/oatk/oatk.asm.mito.gfa"
    output:
        "results/bandage/mito.jpg"
    threads: config["fastp"]["t"]
    log:
        "logs/bandage.log",
    resources:
        mem_mb=config['fastp']['mem_mb'],  # access memory from config
    conda:
        "../envs/bandage.yaml"
    shell:
        """
        Bandage image {input} {output}
        """


rule bandage_pltd:
    input:
        mito_gfa = "results/oatk_pltd/oatk.asm.mito.gfa",
        pltd_gfa = "results/oatk_pltd/oatk.asm.pltd.gfa"
    output:
        mito_out = "results/bandage_pltd/mito.jpg",
        pltd_out = "results/bandage_pltd/pltd.jpg"
    threads: config["fastp"]["t"]
    log:
        "logs/bandage.log",
    resources:
        mem_mb=config['fastp']['mem_mb'],  # access memory from config
    conda:
        "../envs/bandage.yaml"
    shell:
        """
        Bandage image {input.mito_gfa} {output.mito_out}
        Bandage image {input.pltd_gfa} {output.pltd_out}
        """