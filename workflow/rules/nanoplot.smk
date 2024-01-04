# This rule runs nanoplot to qc hifi reads. 

rule run_nanoplot:
    input:
        "results/reads/hifi/hifi.fastq.gz"
    output:
        "results/nanoplot/NanoPlot-report.html"
    log:
        "logs/nanoplot.log"
    conda:
        "../envs/nanoplot.yaml"
    shell:
        """
        NanoPlot -t {config[nanoplot][t]} --fastq {input} --loglength -o results/nanoplot --plots dot --verbose >> {log} 2>&1
        """