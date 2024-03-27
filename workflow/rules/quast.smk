# rule quast for assembly QC

rule quast:
    input:
        completion = checkpoints.yahs_completed.get().completion_marker,
        assemblies = expand("results/assemblies/{file}.fa", file=assembly_files)
    output:
        "results/quast/report.html",
    threads: config["quast"]["t"]
    log:
        "logs/quast.log",
    resources:
        mem_mb=config['quast']['mem_mb'],  # access memory from config
    conda:
        "../envs/quast.yaml"
    params:
        optional_params=" ".join(
            f"{k} {v}" for k, v in config["quast"]["optional_params"].items() if v
        ),
    shell:
        """ 
        quast {input.assemblies} -o results/quast  -t {threads} {params.optional_params}
        """