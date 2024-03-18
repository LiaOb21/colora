# This rule runs the decontamination pipeline fcs-gx on the hifiasm primary assembly
# Github page of the pipeline: https://github.com/ncbi/fcs

rule fcsgx:
    input:
        fasta="results/hifiasm/asm.{hap}.fa",
    params:
        ncbi_tax_id=config["fcsgx"]["ncbi_tax_id"],
        path_to_gx_db=config["fcsgx"]["path_to_gx_db"],
    output:
        clean_fasta="results/ncbi_fcsgx_{hap}/asm_clean.fa",
        action_report="results/ncbi_fcsgx_{hap}/out/asm.{params.ncbi_tax_id}.fcs_gx_report.txt",
        contaminants="results/ncbi_fcsgx_{hap}/contaminants.fa"
    threads: config["hifiasm"]["t"]
    log:
        "logs/fcsgx_{hap}.log",
    resources:
        mem_mb=config['fcsgx']['mem_mb'],  # access memory from config
    conda:
        "../envs/fcsgx.yaml"
    shell:
        """
        run_gx.py --fasta {input.fasta} --out-dir results/ncbi_fcsgx_{wildcards.hap}/out --gx-db {params.path_to_gx_db} --tax-id {params.ncbi_tax_id} --debug >> {log} 2>&1
        cat {input.fasta} | gx clean-genome --action-report {output.action_report} --output {output.clean_fasta} --contam-fasta-out {output.contaminants} >> {log} 2>&1
        """
