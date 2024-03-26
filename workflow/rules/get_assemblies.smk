include: "common.smk"

rule get_assemblies:
    output:
        directory("results/assemblies")
    log:
        "logs/get_assemblies.log",
    resources:
        mem_mb=config['arima']['mem_mb'],  # access memory from config
    conda:
        "../envs/basic.yaml"
    shell:
        """
        mkdir -p {output}
        """