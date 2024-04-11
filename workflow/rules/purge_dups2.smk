# This rule uses purge_dups purge haplotigs and overlaps from the alternate assembly produced by hifiasm


rule purge_dups_alt:
    input:
        reads="results/reads/hifi/hifi.fastq.gz",
        fasta="results/hifiasm/asm.alternate.fa",
        hap_fa_in="results/purge_dups/hap.fa",
    output:
        merged_fasta="results/purge_dups_alt/merged.fa",
        paf="results/purge_dups_alt/hifi_vs_alternate_contigs.paf.gz",
        cutoffs="results/purge_dups_alt/cutoffs",
        pcbstat_stat="results/purge_dups_alt/PB.stat",
        pcbstat_cov="results/purge_dups_alt/PB.base.cov",
        pcbstat_cov_wig="results/purge_dups_alt/PB.cov.wig",
        split_fa="results/purge_dups_alt/asm.alternate.split",
        self_paf="results/purge_dups_alt/asm.alternate.split.self.paf.gz",
        dups_bed="results/purge_dups_alt/dups.bed",
        hap_fa="results/purge_dups_alt/hap.fa",
        purged_fasta="results/purge_dups_alt/asm.alternate_purged.fa",
        hist_plot="results/purge_dups_alt/hist.out.png",
    threads: config["minimap2"]["t"]
    log:
        "logs/purge_dups_alt.log",
    resources:
        mem_mb=config['purge_dups']['mem_mb'],  # access memory from config
    conda:
        "../envs/purge_dups.yaml"
    shell:
        """
        cat {input.fasta} {input.hap_fa_in} > merged.fa 
        echo "results/purge_dups/hap.fa and results/hifiasm/asm.alternate.fa merged" >> {log}
        (minimap2 -xasm20 merged.fa {input.reads} -t {threads} | gzip -c - > hifi_vs_alternate_contigs.paf.gz) 2>> {log}
        pbcstat hifi_vs_alternate_contigs.paf.gz 2>> {log}
        (calcuts PB.stat > cutoffs) 2>> {log}
        (split_fa merged.fa > asm.alternate.split) 2>> {log}
        (minimap2 -xasm5 -DP asm.alternate.split asm.alternate.split -t {threads} | gzip -c - > asm.alternate.split.self.paf.gz) 2>> {log} 
        (purge_dups -2 -T cutoffs -c PB.base.cov asm.alternate.split.self.paf.gz > dups.bed) 2>> {log}
        hist_plot.py -c cutoffs PB.stat hist.out.png 2>> {log}

        mkdir -p results/purge_dups/ 
        mv merged.fa {output.merged_fasta} 
        mv hifi_vs_alternate_contigs.paf.gz {output.paf} 
        mv PB.stat {output.pcbstat_stat} 
        mv PB.base.cov {output.pcbstat_cov} 
        mv PB.cov.wig {output.pcbstat_cov_wig} 
        mv cutoffs {output.cutoffs} 
        mv asm.alternate.split {output.split_fa}
        mv asm.alternate.split.self.paf.gz {output.self_paf}
        mv dups.bed {output.dups_bed}
        mv hap.fa {output.hap_fa}
        mv purged.fa {output.purged_fasta}
        mv hist.out.png {output.hist_plot}
        """
