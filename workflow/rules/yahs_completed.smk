# add a checkpoint to run the QC after all the assemblies have been symlinked to results/assemblies

rule all_yahs_completed:
    input:
        get_yahs_output()
    output:
        "results/yahs_all_done.txt"
    log:
        "logs/all_yahs_completed.log",
    conda:
        "../envs/basic.yaml"
    shell:
        """
        ls {input} > {output}
        """