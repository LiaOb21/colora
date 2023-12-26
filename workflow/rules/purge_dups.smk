# This rule uses purge_dups purge haplotigs and overlaps in the assembly produced by hifiasm

rule run_purge_dups:
    input:
        reads = "results/reads/hifi/hifi.fastq.gz",
        fasta = "results/assemblies/hifiasm/hifiasm.asm.bp.p_ctg.fa"
    output:
        paf = "results/assemblies/purge_dups/hifi_vs_hifiasm_contigs.paf.gz",
        pcbstat_stat = "results/assemblies/purge_dups/PB.stat",
        pcbstat_cov = "results/assemblies/purge_dups/PB.base.cov",
        split_fa = "results/assemblies/purge_dups/hifiasm.asm.split",
        self_paf = "results/assemblies/purge_dups/hifiasm.split.self.paf.gz",
        dups_bed = "results/assemblies/purge_dups/dups.bed",
        log = "results/assemblies/purge_dups/purge_dups.log",
        purged_fasta = "results/assemblies/purge_dups/purged.fa",
        fasta_rename = "results/assemblies/purge_dups/hifiasm.asm.purged.fa"
    log:
        "logs/purge_dups.log"
    conda:
        "../envs/purge_dups.yaml"
    shell:
        """
        minimap2 -xasm20 {input.fasta} {input.reads} -t {config[minimap2][t]} | gzip -c - > {output.paf}
        pbcstat {output.paf} 
        calcuts {output.pcbstat_stat} > cutoffs 2>calcults.log
        split_fa {input.fasta} > {output.split_fa}
        minimap2 -xasm5 -DP {output.split_fa} {output.split_fa} -t {config[minimap2][t]} | gzip -c - > {output.self_paf}
        purge_dups -2 -T cutoffs -c {output.pcbstat_cov} {output.self_paf} > {output.dups_bed} 2> {output.log}
        get_seqs -e {output.dups_bed} {input.fasta} > {output.purged_fasta}

        mkdir -p results/assemblies/purge_dups/
        mv {output.paf} results/assemblies/purge_dups/
        mv PB.* results/assemblies/purge_dups/
        mv calcults.log results/assemblies/purge_dups/
        mv cutoffs results/assemblies/purge_dups/
        mv {output.split_fa} results/assemblies/purge_dups/
        mv {output.self_paf} results/assemblies/purge_dups/
        mv {output.dups_bed} results/assemblies/purge_dups/
        mv {output.log} results/assemblies/purge_dups/
        mv {output.purged_fasta} {output.fasta_rename}
        mv {output.fasta_rename} results/assemblies/purge_dups/
        """