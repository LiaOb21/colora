# rule busco for assembly QC

rule busco:
    input:
        completion = "results/yahs_all_done.txt",
        assemblies = "results/assemblies/{file}.fa"
    output:
        dir = directory("results/busco/{file}.fa_busco"),
    params:
        lineage=config["busco"]["lineage"],
    threads: config["busco"]["t"]
    log:
        "logs/busco{file}.log"
    resources:
        mem_mb=config['busco']['mem_mb'],  # access memory from config
    conda:
        "../envs/busco.yaml"
    shell:
        """
        # run busco on each assembly
		busco -i {input.assemblies} -o {output.dir} -m genome -l {params.lineage} -f -c {threads}
        """