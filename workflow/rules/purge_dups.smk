# This rule uses purge_dups purge haplotigs and overlaps from the primary assembly produced by hifiasm


# include common.smk to use get_purge_dups_inputs
include: "common.smk"


rule purge_dups:
    input:
        reads="results/reads/hifi/hifi.fastq.gz",
        fasta=get_purge_dups_inputs(),
    output:
        paf="results/purge_dups/hifi_vs_primary_contigs.paf.gz",
        calcuts_log="results/purge_dups/calcuts.log",
        cutoffs="results/purge_dups/cutoffs",
        pcbstat_stat="results/purge_dups/PB.stat",
        pcbstat_cov="results/purge_dups/PB.base.cov",
        pcbstat_cov_wig="results/purge_dups/PB.cov.wig",
        split_fa="results/purge_dups/asm.primary.split",
        self_paf="results/purge_dups/asm.primary.split.self.paf.gz",
        dups_bed="results/purge_dups/dups.bed",
        log="results/purge_dups/purge_dups.log",
        hap_fa="results/purge_dups/hap.fa",
        purged_fasta="results/purge_dups/asm.primary_purged.fa",
        hist_plot="results/purge_dups/hist.out.png",
    threads: config["minimap2"]["t"]
    log:
        "logs/purge_dups.log",
    resources:
        mem_mb=config['purge_dups']['mem_mb'],  # access memory from config
    conda:
        "../envs/purge_dups.yaml"
    shell:
        """
        minimap2 -xasm20 {input.fasta} {input.reads} -t {threads} | gzip -c - > hifi_vs_primary_contigs.paf.gz 
        pbcstat hifi_vs_primary_contigs.paf.gz 
        calcuts PB.stat > cutoffs 2>calcults.log 
        split_fa results/hifiasm/asm.primary.fa > asm.primary.split 
        minimap2 -xasm5 -DP asm.primary.split asm.primary.split -t {threads} | gzip -c - > asm.primary.split.self.paf.gz 
        purge_dups -2 -T cutoffs -c PB.base.cov asm.primary.split.self.paf.gz > dups.bed 2> purge_dups.log 
        get_seqs -e dups.bed results/hifiasm/asm.primary.fa > purged.fa
        hist_plot.py -c cutoffs PB.stat hist.out.png 

        mkdir -p results/purge_dups/ 
        mv hifi_vs_primary_contigs.paf.gz {output.paf}
        mv PB.stat {output.pcbstat_stat}
        mv PB.base.cov {output.pcbstat_cov}
        mv PB.cov.wig {output.pcbstat_cov_wig}
        mv calcults.log {output.calcuts_log}
        mv cutoffs {output.cutoffs}
        mv asm.primary.split {output.split_fa}
        mv asm.primary.split.self.paf.gz {output.self_paf}
        mv dups.bed {output.dups_bed}
        mv purge_dups.log {output.log}
        mv hap.fa {output.hap_fa}
        mv purged.fa {output.purged_fasta}
        mv hist.out.png {output.hist_plot}
        """
