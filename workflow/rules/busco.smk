# rule busco for assembly QC

rule busco:
    input:
        unpack(get_assemblies_QC)
    output:
        dir = directory("results/busco"),
    params:
        lineage=config["busco"]["lineage"],
        labels = list(get_assemblies_QC().keys())
    threads: config["busco"]["t"]
    log:
        "logs/busco.log"
    resources:
        mem_mb=config['busco']['mem_mb'],  # access memory from config
    conda:
        "../envs/busco.yaml"
    shell:
        """
        mkdir -p {output}
        labels=({params.labels})
        count=0
        for v in {input} ; do
            k=${{labels[$count]}}
            /usr/bin/time -v busco -i $v -o {output}/$k -m genome -l {params.lineage} -f -c {threads} >> {log} 2>&1
            count=$(( $count + 1 ))
        done
        """