# This rule runs nanoplot to qc hifi reads. 

rule nanoplot:
    input:
        "results/reads/hifi/hifi.fastq.gz"
    output:
        "results/nanoplot/NanoPlot-report.html"
    threads: config['nanoplot']['t']
    log:
        "logs/nanoplot.log"
    conda:
        "../envs/nanoplot.yaml"
    shell:
        """
        NanoPlot -t {threads} --fastq {input} --loglength -o results/nanoplot --plots dot --verbose >> {log} 2>&1
        """