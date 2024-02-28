# This rule runs the decontamination pipeline fcs-gx on the hifiasm primary assembly
# Github page of the pipeline: https://github.com/ncbi/fcs


rule fcsgx:
    input:
        fasta = "results/hifiasm/hifiasm.asm.p_ctg.fa"
    params:
        ncbi_tax_id = config["fcsgx"]["ncbi_tax_id"],
        action_report = f"results/ncbi_fcsgx/hifiasm.asm.p_ctg.fa.{config['fcsgx']['ncbi_tax_id']}.fcs_gx_report.txt",
        path_to_gx_db = config["fcsgx"]["path_to_gx_db"],
        contaminants_output_name = config["fcsgx"]["contaminants_output_name"]
    output:
        action_report = "{params.action_report}",
        clean_fasta="results/ncbi_fcsgx/hifiasm.asm.p_ctg_clean.fa",
    threads: config['hifiasm']['t']
    log:
        "logs/hifiasm.log"
    conda:
        "../envs/fcsgx.yaml"
    shell:
        """
        run_gx.py --fasta {input.fasta} --out-dir results/ncbi_fcsgx --gx-db {params.path_to_gx_db} --tax-id {params.ncbi_tax_id} --debug >> {log} 2>&1
        cat {input.fasta} | python3 scripts/fcs.py clean genome --action-report {params.action_report} --output {output.clean_fasta} --contam-fasta-out {params.contaminants_output_name} >> {log} 2>&1
        """