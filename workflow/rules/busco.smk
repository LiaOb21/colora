# rule busco for assembly QC

rule busco:
    input:
        completion = checkpoints.yahs_completed.get().completion_marker,
        assemblies = "results/assemblies/{file}.fa"
    output:
        dir = directory("results/busco/{file}.fa_busco"),
    params:
        lineage=config["busco"]["lineage"],
    threads: config["busco"]["t"]
    log:
        "logs/busco_{file}.log"
    resources:
        mem_mb=config['busco']['mem_mb'],  # access memory from config
    conda:
        "../envs/busco.yaml"
    shell:
        """
        # run busco on each assembly
        #assembly_name=$(basename {input.assemblies})
		busco -i {input.assemblies} -o {output.dir}/{wildcards.file}_busco -m genome -l {params.lineage} -f -c {threads}
        """