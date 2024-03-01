# This rule pre-processes hic reads using fastp.

import glob


rule fastp:
    input:
        forward_in=expand(
            "{hic_path}{sample}_1.fastq.gz",
            hic_path=config["hic_path"],
            sample=glob_wildcards(config["hic_path"] + "{sample}_1.fastq.gz").sample,
        ),
        reverse_in=expand(
            "{hic_path}{sample}_2.fastq.gz",
            hic_path=config["hic_path"],
            sample=glob_wildcards(config["hic_path"] + "{sample}_2.fastq.gz").sample,
        ),
    output:
        forward_out="results/fastp/{sample}_trim_1.fastq.gz",
        reverse_out="results/fastp/{sample}_trim_2.fastq.gz",
        json="results/fastp/{sample}_report_fastp.HiC.json",
        html="results/fastp/{sample}_report_fastp.HiC.html",
    params:
        cut_window_size=config["fastp"]["cut_window_size"],
        cut_mean_quality=config["fastp"]["cut_mean_quality"],
        optional_params=" ".join(
            f"{k} {v}" for k, v in config["fastp"]["optional_params"].items() if v
        ),
    threads: config["fastp"]["t"]
    log:
        "logs/{sample}_fastp.log",
    conda:
        "../envs/fastp.yaml"
    shell:
        """
        fastp -p -i {input.forward_in} -I {input.reverse_in} -o {output.forward_out} -O {output.reverse_out} --cut_window_size {params.cut_window_size} --cut_mean_quality {params.cut_mean_quality} --json {output.json} --html {output.html} --thread {threads} {params.optional_params} >> {log} 2>&1
        """
