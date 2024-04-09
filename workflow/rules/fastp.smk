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
    threads: config["fastp"]["t"]
    log:
        "logs/fastp.log",
    resources:
        mem_mb=config['fastp']['mem_mb'],  # access memory from config
    conda:
        "../envs/fastp.yaml"
    shell:
        """
        fastp -p -i {input.forward_in} -I {input.reverse_in} -o {output.forward_out} -O {output.reverse_out} --cut_front --cut_front_window_size 5  --detect_adapter_for_pe --json {output.json} --html {output.html} --thread {threads} >> {log} 2>&1
        """
