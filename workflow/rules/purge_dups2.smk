# This rule uses purge_dups purge haplotigs and overlaps from the alternate assembly produced by hifiasm


rule purge_dups_alt:
    input:
        reads="results/reads/hifi/hifi.fastq.gz",
        fasta="results/hifiasm/asm.alternate.fa",
        hap_fa_in="results/purge_dups/hap.fa",
        asm_dir = "results/assemblies"
    output:
        merged_fasta="results/purge_dups_alt/merged.fa",
        paf="results/purge_dups_alt/hifi_vs_alternate_contigs.paf.gz",
        calcuts_log="results/purge_dups_alt/calcuts.log",
        cutoffs="results/purge_dups_alt/cutoffs",
        pcbstat_stat="results/purge_dups_alt/PB.stat",
        pcbstat_cov="results/purge_dups_alt/PB.base.cov",
        pcbstat_cov_wig="results/purge_dups_alt/PB.cov.wig",
        split_fa="results/purge_dups_alt/asm.alternate.split",
        self_paf="results/purge_dups_alt/asm.alternate.split.self.paf.gz",
        dups_bed="results/purge_dups_alt/dups.bed",
        log="results/purge_dups_alt/purge_dups.log",
        hap_fa="results/purge_dups_alt/hap.fa",
        purged_fasta="results/purge_dups_alt/asm.alternate_purged.fa",
        hist_plot="results/purge_dups_alt/hist.out.png",
        link="results/assemblies/purged_alternate.fa"
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
        minimap2 -xasm20 merged.fa {input.reads} -t {threads} | gzip -c - > hifi_vs_alternate_contigs.paf.gz 
        pbcstat hifi_vs_alternate_contigs.paf.gz 
        calcuts PB.stat > cutoffs 2>calcults.log
        split_fa merged.fa > asm.alternate.split 
        minimap2 -xasm5 -DP asm.alternate.split asm.alternate.split -t {threads} | gzip -c - > asm.alternate.split.self.paf.gz 
        purge_dups -2 -T cutoffs -c PB.base.cov asm.alternate.split.self.paf.gz > dups.bed 2> purge_dups.log
        get_seqs -e dups.bed merged.fa > purged.fa
        hist_plot.py -c cutoffs PB.stat hist.out.png 

        mkdir -p results/purge_dups/ 
        mv merged.fa {output.merged_fasta} 
        mv hifi_vs_alternate_contigs.paf.gz {output.paf} 
        mv PB.stat {output.pcbstat_stat} 
        mv PB.base.cov {output.pcbstat_cov} 
        mv PB.cov.wig {output.pcbstat_cov_wig} 
        mv calcults.log {output.calcuts_log} 
        mv cutoffs {output.cutoffs} 
        mv asm.alternate.split {output.split_fa}
        mv asm.alternate.split.self.paf.gz {output.self_paf}
        mv dups.bed {output.dups_bed}
        mv purge_dups.log {output.log}
        mv hap.fa {output.hap_fa}
        mv purged.fa {output.purged_fasta}
        mv hist.out.png {output.hist_plot}

        ln -srn {output.purged_fasta} {output.link}
        """
