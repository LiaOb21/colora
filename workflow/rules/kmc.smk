# This rule generates kmers from hifi raw data

rule run_kmc:
    input:
        "results/reads/hifi/hifi.fastq.gz"
    output:
        out_kmc_pre = "results/kmc/output.kmc_pre",
        out_kmc_suf = "results/kmc/output.kmc_suf",
        hist = "results/kmc/out.hist"
    log:
        "logs/kmc.log"
    conda:
        "../envs/kmc.yaml"
    shell:
        """
        mkdir -p results/kmc/temp
        kmc -k{config[kmc][k]} -t{config[kmc][t]} -ci{config[kmc][ci]} -cs{config[kmc][cs]} -fq {input} results/kmc/output results/kmc/temp/
        kmc_tools transform results/kmc/output histogram {output.hist} -cx{config[kmc_tools][cx]}
        """