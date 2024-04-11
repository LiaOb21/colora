# This rule uses purge_dups purge haplotigs and overlaps from the primary assembly produced by hifiasm

rule purge_dups:
    input:
        reads = "results/reads/hifi/hifi.fastq.gz",
        fasta = get_purge_dups_inputs(),
    output:
        paf="results/purge_dups/hifi_vs_primary_contigs.paf.gz",
        cutoffs="results/purge_dups/cutoffs",
        pcbstat_stat="results/purge_dups/PB.stat",
        pcbstat_cov="results/purge_dups/PB.base.cov",
        pcbstat_cov_wig="results/purge_dups/PB.cov.wig",
        split_fa="results/purge_dups/asm.primary.split",
        self_paf="results/purge_dups/asm.primary.split.self.paf.gz",
        dups_bed="results/purge_dups/dups.bed",
        hap_fa="results/purge_dups/hap.fa",
        purged_fasta="results/purge_dups/asm.primary_purged.fa",
        hist_plot="results/purge_dups/hist.out.png",
        link = "results/assemblies/purged_primary.fa",
    threads: config["minimap2"]["t"]
    log:
        "logs/purge_dups.log",
    resources:
        mem_mb=config['purge_dups']['mem_mb'],  # access memory from config
    conda:
        "../envs/purge_dups.yaml"
    shell:
        """
        (minimap2 -xasm20 {input.fasta} {input.reads} -t {threads} | gzip -c - > hifi_vs_primary_contigs.paf.gz) 2>> {log}
        pbcstat hifi_vs_primary_contigs.paf.gz 2>> {log}
        (calcuts PB.stat > cutoffs) 2>> {log} 
        (split_fa results/hifiasm/asm.primary.fa > asm.primary.split) 2>> {log} 
        (minimap2 -xasm5 -DP asm.primary.split asm.primary.split -t {threads} | gzip -c - > asm.primary.split.self.paf.gz) 2>> {log}
        (purge_dups -2 -T cutoffs -c PB.base.cov asm.primary.split.self.paf.gz > dups.bed) 2>> {log} 
        (get_seqs -e dups.bed results/hifiasm/asm.primary.fa > purged.fa) 2>> {log}
        hist_plot.py -c cutoffs PB.stat hist.out.png 2>> {log}

        mkdir -p results/purge_dups/ 
        mv hifi_vs_primary_contigs.paf.gz {output.paf}
        mv PB.stat {output.pcbstat_stat}
        mv PB.base.cov {output.pcbstat_cov}
        mv PB.cov.wig {output.pcbstat_cov_wig}
        mv cutoffs {output.cutoffs}
        mv asm.primary.split {output.split_fa}
        mv asm.primary.split.self.paf.gz {output.self_paf}
        mv dups.bed {output.dups_bed}
        mv hap.fa {output.hap_fa}
        mv purged.fa {output.purged_fasta}
        mv hist.out.png {output.hist_plot}

        ln -srn {output.purged_fasta} {output.link}
        """
