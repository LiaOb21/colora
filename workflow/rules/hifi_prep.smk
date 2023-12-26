import glob

# This rule joins together HiFi reads in the case in which 
# we have different files from different runs. 
# If the file is only one this rule is skipped

rule hifi_prep:
    input:
        files = expand("{hifi_path}{file}.fastq.gz", hifi_path=config["hifi_path"], file=glob_wildcards(config["hifi_path"] + "{file}.fastq.gz").file)
    output:
        hifi = "results/reads/hifi/hifi.fastq.gz"
    params:
        num_files = len(glob.glob(config["hifi_path"] + "*.fastq.gz"))
    log:
        "logs/hifi_prep.log"
    shell:
        """
        if [ {params.num_files} -gt 1 ]; then
            zcat {input.files} | gzip > {output.hifi}
        else
            cp {input.files[0]} {output.hifi}
        fi
        """