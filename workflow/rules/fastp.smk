# This rule pre-processes hic reads using fastp.

import glob


rule fastp:
    input:
        forward_in=f"{config["hic_path"]}{sample}_1.fastq.gz",
        reverse_in=f"{config["hic_path"]}{sample}_2.fastq.gz",
    output:
        forward_out="results/fastp/hic_trim_1.fastq.gz",
        reverse_out="results/fastp/hic_trim_2.fastq.gz",
        json="results/fastp/hic_report_fastp.HiC.json",
        html="results/fastp/hic_report_fastp.HiC.html",
    params:
        optional_params=" ".join(
            f"{k} {v}" for k, v in config["fastp"]["optional_params"].items() if v
        ),    
    threads: config["medium"]["t"]
    log:
        "logs/fastp.log",
    resources:
        mem_mb=config["medium"]["mem_mb"],  # access memory from config
    conda:
        "../envs/fastp.yaml"
    shell:
        """
        /usr/bin/time -v fastp -p -i {input.forward_in} -I {input.reverse_in} -o {output.forward_out} -O {output.reverse_out} --detect_adapter_for_pe --json {output.json} --html {output.html} --thread {threads} {params.optional_params} >> {log} 2>&1
        """
