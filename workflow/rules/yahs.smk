# This rule uses yahs to scaffold the contigs from hifiasm_p_purged.fa

import glob

rule yahs:
    input:
        REF="results/bwa_index_{hap}/asm.fa",
        bam="results/arima_mapping_pipeline_{hap}/REP_DIR/{sample}_rep1.bam",
        asm_dir = "results/assemblies"
    output:
        scaffolds="results/yahs_{hap}/asm_yahs_{sample}_scaffolds_final.fa",
        link="results/assemblies/yahs_{hap}_{sample}.fa",
    log:
        "logs/yahs_{hap}_{sample}.log",
    resources:
        mem_mb=config['arima']['mem_mb'],  # access memory from config
    conda:
        "../envs/yahs.yaml"
    params:
        optional_params=" ".join(
            f"{k} {v}" for k, v in config["yahs"]["optional_params"].items() if v
        ),
    shell:
        """ 
        yahs {input.REF} {input.bam} -o results/yahs_{wildcards.hap}/asm_yahs_{wildcards.sample} {params.optional_params} >> {log} 2>&1
        ln -srn {output.scaffolds} {output.link}
        """
