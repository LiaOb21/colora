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
    threads: config["high"]["t"]
    log:
        "logs/purge_dups.log",
    resources:
        mem_mb=config["high"]["mem_mb"],  # access memory from config
    conda:
        "../envs/purge_dups.yaml"
    shell:
        """
        (minimap2 -xasm20 {input.fasta} {input.reads} -t {threads} | gzip -c - > hifi_vs_primary_contigs.paf.gz) >> {log} 2>&1
        pbcstat hifi_vs_primary_contigs.paf.gz >> {log} 2>&1
        (calcuts PB.stat > cutoffs) >> {log} 2>&1 
        (split_fa results/hifiasm/asm.primary.fa > asm.primary.split) >> {log} 2>&1
        (minimap2 -xasm5 -DP asm.primary.split asm.primary.split -t {threads} | gzip -c - > asm.primary.split.self.paf.gz) >> {log} 2>&1
        (purge_dups -2 -T cutoffs -c PB.base.cov asm.primary.split.self.paf.gz > dups.bed) >> {log} 2>&1
        (get_seqs -e dups.bed results/hifiasm/asm.primary.fa > purged.fa) >> {log} 2>&1
        hist_plot.py -c cutoffs PB.stat hist.out.png >> {log} 2>&1

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
