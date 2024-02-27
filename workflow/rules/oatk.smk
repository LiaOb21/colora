# This rule runs oatk to extract organelles reads from the hifi reads. 

rule run_oatk:
    input:
        "results/reads/hifi/hifi.fastq.gz"
    output:
        mito = "results/oatk/oatk.asm.mito.ctg.fasta"
    log:
        "logs/oatk.log"
    conda:
        "../envs/oatk.yaml"
    threads: config['oatk']['t']
    params:
        k = config['oatk']['k'],
        c = config['oatk']['c'],
        m = config['oatk']['m'],
        optional_params=" ".join(
            f"{k} {v}" for k, v in config["oatk"]["optional_params"].items() if v
        )
    shell:
        """
        oatk -k {params.k} -c {params.c} -t {threads} -m {params.m} {params.optional_params} {input} >> {log} 2>&1
        mv oatk.asm* results/oatk/ >> {log} 2>&1
        """