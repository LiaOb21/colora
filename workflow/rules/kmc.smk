# This rule generates kmers from hifi raw data


rule kmc:
    input:
        "results/reads/hifi/hifi.fastq.gz",
    output:
        out_kmc_pre="results/kmc/output.kmc_pre",
        out_kmc_suf="results/kmc/output.kmc_suf",
        hist="results/kmc/out.hist",
    threads: config["kmc"]["t"]
    params:
        k=config["kmc"]["k"],
        ci=config["kmc"]["ci"],
        cs=config["kmc"]["cs"],
        cx=config["kmc_tools"]["cx"],
    log:
        "logs/kmc.log",
    resources:
        mem_mb=config['kmc']['mem_mb'],  # access memory from config
    conda:
        "../envs/kmc.yaml"
    shell:
        """
        mkdir -p results/kmc/temp >> {log} 2>&1
        /usr/bin/time -v kmc -k{params.k} -t{threads} -ci{params.ci} -cs{params.cs} -fq {input} results/kmc/output results/kmc/temp/ >> {log} 2>&1
        /usr/bin/time -v kmc_tools transform results/kmc/output histogram {output.hist} -cx{params.cx} >> {log} 2>&1
        """
