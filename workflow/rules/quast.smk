# rule quast for assembly QC

rule quast:
    input:
        unpack(get_assemblies_QC)
    output:
        directory("results/quast"),
    threads: config["medium"]["t"]
    log:
        "logs/quast.log",
    resources:
        mem_mb=config["medium"]["mem_mb"],  # access memory from config
    conda:
        "../envs/quast.yaml"
    params:
        optional_params=" ".join(
            f"{k} {v}" for k, v in config["quast"]["optional_params"].items() if v
        ),
        labels = ",".join(get_assemblies_QC().keys())
    shell:
        """
        quast {input} -l {params.labels} -o {output}  -t {threads} {params.optional_params} >> {log} 2>&1
        """
