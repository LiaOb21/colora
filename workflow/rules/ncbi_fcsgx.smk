# This rule runs the decontamination pipeline fcs-gx on the hifiasm primary assembly
# Github page of the pipeline: https://github.com/ncbi/fcs


rule fcsgx:
    input:
        fasta="results/hifiasm/hifiasm.asm.p_ctg.fa",
    params:
        ncbi_tax_id=config["fcsgx"]["ncbi_tax_id"],
        action_report_name=config["fcsgx"]["action_report_name"],
        path_to_gx_db=config["fcsgx"]["path_to_gx_db"],
        contaminants_output_name=config["fcsgx"]["contaminants_output_name"],
    output:
        clean_fasta="results/ncbi_fcsgx/hifiasm.asm.p_ctg_clean.fa",
    threads: config["hifiasm"]["t"]
    log:
        "logs/fcsgx.log",
    resources:
        mem_mb=config['fcsgx']['mem_mb'],  # access memory from config
    conda:
        "../envs/fcsgx.yaml"
    shell:
        """
        run_gx.py --fasta {input.fasta} --out-dir results/ncbi_fcsgx/out --gx-db {params.path_to_gx_db} --tax-id {params.ncbi_tax_id} --debug >> {log} 2>&1
        cat {input.fasta} | gx clean-genome --action-report results/ncbi_fcsgx/out/{params.action_report_name} --output {output.clean_fasta} --contam-fasta-out results/ncbi_fcsgx/{params.contaminants_output_name} >> {log} 2>&1
        """
