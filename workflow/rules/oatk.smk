# This rule runs oatk to extract organelles reads from the hifi reads. 

rule run_oatk:
    input:
        "results/reads/hifi/hifi.fastq.gz"
    output:
        mito = "results/oatk/oatk.asm.mito.ctg.fasta",
        pltd = "results/oatk/oatk.asm.pltd.ctg.fasta"
    conda:
        "../envs/oatk.yaml"
    shell:
        """
        oatk -k {config[oatk][k]} -c {config[oatk][c]} -t {config[oatk][t]} -m {config[oatk][m]} -p {config[oatk][p]} {input}
        mv oatk.asm* results/oatk/
        """