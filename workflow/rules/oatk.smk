# This rule runs oatk to extract organelles reads from the hifi reads. 

rule run_oatk:
    input:
        "results/reads/hifi/hifi.fastq.gz"
    output:
        mito = "results/oatk/oatk.asm.mito.ctg.fasta",
        pltd = "results/oatk/oatk.asm.pltd.ctg.fasta"
    log:
        "logs/oatk.log"
    conda:
        "../envs/oatk.yaml"
    params:
        optional_params=" ".join(
            f"{k} {v}" for k, v in config["oatk"]["optional_params"].items() if v
        )
    shell:
        """
        oatk -k {config[oatk][k]} -c {config[oatk][c]} -t {config[oatk][t]} -m {config[oatk][m]} {params.optional_params} {input} >> {log} 2>&1
        mv oatk.asm* results/oatk/ >> {log} 2>&1
        """