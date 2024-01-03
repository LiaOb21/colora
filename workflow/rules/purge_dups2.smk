# This rule uses purge_dups purge haplotigs and overlaps from the alternate assembly produced by hifiasm

rule run_purge_dups_alt:
    input:
        reads = "results/reads/hifi/hifi.fastq.gz",
        fasta = "results/hifiasm/hifiasm.asm.a_ctg.fa",
        hap_fa_in = "results/purge_dups/hap.fa"
    output:
        merged_fasta = "results/purge_dups_alt/merged.fa",
        paf = "results/purge_dups_alt/hifi_vs_hifiasm_contigs.paf.gz",
        calcuts_log = "results/purge_dups_alt/calcuts.log",
        cutoffs = "results/purge_dups_alt/cutoffs",
        pcbstat_stat = "results/purge_dups_alt/PB.stat",
        pcbstat_cov = "results/purge_dups_alt/PB.base.cov",
        pcbstat_cov_wig = "results/purge_dups_alt/PB.cov.wig",
        split_fa = "results/purge_dups_alt/hifiasm.asm.split",
        self_paf = "results/purge_dups_alt/hifiasm.split.self.paf.gz",
        dups_bed = "results/purge_dups_alt/dups.bed",
        log = "results/purge_dups_alt/purge_dups.log",
        hap_fa = "results/purge_dups_alt/hap.fa",
        purged_fasta = "results/purge_dups_alt/hifiasm_a_purged.fa",
        hist_plot = "results/purge_dups_alt/hist.out.png"
    conda:
        "../envs/purge_dups.yaml"
    log:
        "logs/purge_dups_alt.log"
    shell:
        """
        cat {input.fasta} {input.hap_fa_in} > merged.fa >> {log} 2>&1
        minimap2 -xasm20 merged.fa {input.reads} -t {config[minimap2][t]} | gzip -c - > hifi_vs_hifiasm_contigs.paf.gz >> {log} 2>&1
        pbcstat hifi_vs_hifiasm_contigs.paf.gz >> {log} 2>&1
        calcuts PB.stat > cutoffs 2>calcults.log
        split_fa merged.fa > hifiasm.asm.split >> {log} 2>&1
        minimap2 -xasm5 -DP hifiasm.asm.split hifiasm.asm.split -t {config[minimap2][t]} | gzip -c - > hifiasm.split.self.paf.gz >> {log} 2>&1
        purge_dups -2 -T cutoffs -c PB.base.cov hifiasm.split.self.paf.gz > dups.bed 2> purge_dups.log
        get_seqs -e dups.bed results/hifiasm/hifiasm.asm.a_ctg.fa > hifiasm_a_purged.fa >> {log} 2>&1
        hist_plot.py -c cutoffs PB.stat hist.out.png >> {log} 2>&1

        mkdir -p results/purge_dups/ >> {log} 2>&1
        mv merged.fa {output.merged_fasta} >> {log} 2>&1
        mv hifi_vs_hifiasm_contigs.paf.gz {output.paf} >> {log} 2>&1
        mv PB.stat {output.pcbstat_stat} >> {log} 2>&1
        mv PB.base.cov {output.pcbstat_cov} >> {log} 2>&1
        mv PB.cov.wig {output.pcbstat_cov_wig} >> {log} 2>&1
        mv calcults.log {output.calcuts_log} >> {log} 2>&1
        mv cutoffs {output.cutoffs} >> {log} 2>&1
        mv hifiasm.asm.split {output.split_fa} >> {log} 2>&1
        mv hifiasm.split.self.paf.gz {output.self_paf} >> {log} 2>&1
        mv dups.bed {output.dups_bed} >> {log} 2>&1
        mv purge_dups.log {output.log} >> {log} 2>&1
        mv hap.fa {output.hap_fa} >> {log} 2>&1
        mv purged.fa {output.purged_fasta} >> {log} 2>&1
        mv hist.out.png {output.hist_plot} >> {log} 2>&1
        """