# This rule runs nanoplot to qc hifi reads. 

rule run_nanoplot:
    input:
        "results/reads/hifi/hifi.fastq.gz"
    output:
        "results/nanoplot"
    log:
        "logs/nanoplot.log"
    conda:
        "../envs/nanoplot.yaml"
    shell:
        """
        NanoPlot -t {config[nanoplot][t]} --fastq {input} --loglength -o {output} --plots dot --verbose
        """