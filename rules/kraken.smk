rule getKrakenTaxo:
    output:
        directory("data/KrakenDB/")
    conda:
        "../envs/kraken2.yaml"
    threads:
        16
    log:
        "logs/Kraken2/getKrakenTaxo.log"
    shell:
        "kraken2-build "
        "--threads {threads} "
        "--download-taxonomy --use-ftp "
        "--db {output} > {log}"


rule getKrakenDB:
    input:
        "data/KrakenDB/"
    output:
        directory("data/KrakenDB/library/{species}")
    conda:
        "../envs/kraken2.yaml"
    threads:
        16
    log:
        "logs/Kraken2/getKrakenDB_{species}.log"
    shell:
        "kraken2-build "
        "--threads {threads} "
        "--download-library {wildcards.species} --use-ftp "
        "--db {input[0]} >> {log}"


rule buildKraken:
    input:
        expand(
            "data/KrakenDB/library/{species}",
            species=config['data']['kraken_species']
        ),
        db = "data/KrakenDB/"
    output:
        "data/KrakenDB.tkn"
    conda:
        "../envs/kraken2.yaml"
    threads:
        32
    log:
        "logs/Kraken2/buildKraken.log"
    shell:
        "kraken2-build "
        "--threads {threads} "
        "--build "
        "--db {input.db} >> {log} "
        "&& touch {output}"


rule runningKraken:
    input:
        tok = "data/KrakenDB.tkn",
        fq1 = "data/trimmed/{id}_1.fastq.gz",
        fq2 = "data/trimmed/{id}_2.fastq.gz",
        db = "data/KrakenDB/"
    output:
        report = "results/Kraken2/{id}/results.report",
        results = "results/Kraken2/{id}/results.smt",
        cseq_reads1 = "results/Kraken2/{id}/cseq_1.fastq",
        cseq_reads2 = "results/Kraken2/{id}/cseq_2.fastq"
    params:
        creads = "results/Kraken2/{id}/cseq#.fastq"
    conda:
        "../envs/kraken2.yaml"
    threads:
        32
    log:
        "logs/Kraken2/runningKraken_{id}.log"
    shell:
        "kraken2 "
        "--threads {threads} "
        "--db {input.db} "
        "--paired --check-names "
        "--report {output.report} "
        "--use-names "
        "--fastq-input "
        "--gzip-compressed "
        "--classified-out {params.creads} "
        "--output {output.results} "
        "{input.fq1} {input.fq2} "
        "> {log}"
