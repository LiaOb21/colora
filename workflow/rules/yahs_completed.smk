# add a checkpoint to run the QC after all the assemblies have been symlinked to results/assemblies

checkpoint symlinks:
    input:
        expand("results/assemblies/yahs_{hap}_{sample}.fa", sample=samples, hap=hap)
    output:
        completion_marker = "results/yahs_all_done.txt"
    log:
        "logs/all_yahs_completed.log",
    conda:
        "../envs/basic.yaml"
    shell:
        """
        ls {input} > {output.completion_marker}
        """