# add a checkpoint to run the QC after all the assemblies have been symlinked to results/assemblies

checkpoint yahs_done:
    input:
        get_yahs_output()
    output:
        completion_marker = "results/yahs_completed.txt"
    log:
        "logs/yahs_completed.log",
    conda:
        "../envs/basic.yaml"
    shell:
        """
        ls {input} > {output.completion_marker}
        """