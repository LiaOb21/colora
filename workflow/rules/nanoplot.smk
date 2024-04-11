# This rule runs nanoplot to qc hifi reads.


rule nanoplot:
    input:
        "results/reads/hifi/hifi.fastq.gz",
    output:
        "results/nanoplot/NanoPlot-report.html",
    threads: config["medium"]["t"]
    log:
        "logs/nanoplot.log",
    resources:
        mem_mb=config["medium"]["mem_mb"],  # access memory from config
    conda:
        "../envs/nanoplot.yaml"
    shell:
        """
        /usr/bin/time -v NanoPlot -t {threads} --fastq {input} --loglength -o results/nanoplot --plots dot --verbose >> {log} 2>&1
        """
