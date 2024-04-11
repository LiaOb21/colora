# This rule uses yahs to scaffold the contigs from hifiasm_p_purged.fa

import glob

rule yahs:
    input:
        REF="results/bwa_index_{hap}/asm.fa",
        bam="results/arima_mapping_pipeline_{hap}/REP_DIR/paired_mark_dups_final.bam",
    output:
        scaffolds="results/yahs_{hap}/asm_yahs_scaffolds_final.fa",
        link="results/assemblies/yahs_{hap}.fa"
    log:
        "logs/yahs_{hap}.log",
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
        /usr/bin/time -v yahs {input.REF} {input.bam} -o results/yahs_{wildcards.hap}/asm_yahs {params.optional_params} >> {log} 2>&1
        ln -srn {output.scaffolds} {output.link}
        """
