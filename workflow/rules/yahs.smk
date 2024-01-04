# This rule uses yahs to scaffold the contigs from hifiasm_p_purged.fa

import glob

rule run_yahs:
    input:
        REF = "results/purge_dups/hifiasm_p_purged.fa",
        bam="results/arima_mapping_pipeline/REP_DIR/{sample}_rep1.bam"
        #bam=expand("results/arima_mapping_pipeline/REP_DIR/{sample}_rep1.bam", sample=glob_wildcards("results/arima_mapping_pipeline/REP_DIR/{sample}_rep1.bam").sample)
    output:
        scaffolds="results/yahs/hifiasm_p_purged_yahs_{sample}_scaffolds_final.fa"
    log:
        "logs/yahs_{sample}.log"
    conda:
        "../envs/yahs.yaml"
    params:
        optional_params=" ".join(
            f"{k} {v}" for k, v in config["yahs"]["optional_params"].items() if v
        )
    shell:
        """
        samtools faidx {input.REF} 
        mkdir -p results/yahs 
        yahs {input.REF} {input.bam} -o results/yahs/hifiasm_p_purged_yahs {params.optional_params} >> {log} 2>&1
        """