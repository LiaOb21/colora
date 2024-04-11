rule bandage:
    input:
        "results/oatk/oatk.asm.mito.gfa"
    output:
        "results/bandage/mito.jpg"
    threads: config["low"]["t"]
    log:
        "logs/bandage.log",
    resources:
        mem_mb=config["low"]["mem_mb"],  # access memory from config
    conda:
        "../envs/bandage.yaml"
    shell:
        """
        /usr/bin/time -v Bandage image {input} {output} &>> {log}
        """


rule bandage_pltd:
    input:
        mito_gfa = "results/oatk_pltd/oatk.asm.mito.gfa",
        pltd_gfa = "results/oatk_pltd/oatk.asm.pltd.gfa"
    output:
        mito_out = "results/bandage_pltd/mito.jpg",
        pltd_out = "results/bandage_pltd/pltd.jpg"
    threads: config["low"]["t"]
    log:
        "logs/bandage.log",
    resources:
        mem_mb=config["low"]["mem_mb"],   # access memory from config
    conda:
        "../envs/bandage.yaml"
    shell:
        """
        /usr/bin/time -v Bandage image {input.mito_gfa} {output.mito_out} >> {log} 2>&1
        /usr/bin/time -v Bandage image {input.pltd_gfa} {output.pltd_out} >> {log} 2>&1
        """