rule DESeq2_genes:
    """ Differential expression for the different conditions """
    input:
        counts = 'results/kallisto/est_counts.tsv',
        samples = "data/design.tsv"
    output:
        results = directory("results/DESeq2/genes_test"),
    log:
        "logs/DESeq2/genes.log"
    conda:
        "../envs/R.yaml"
    script:
        "../scripts/DESeq2_genes.R"

rule DESeq2_transcripts:
    """ Differential expression for the different conditions """
    input:
        counts = expand("results/kallisto/{id}/abundance.h5",
                            id=simple_id),
        samples = "data/design.tsv"
    output:
        results = directory("results/DESeq2/transcripts")
    params:
        names = expand('{id}', id=simple_id)
    log:
        "logs/DESeq2/transcripts.log"
    conda:
        "../envs/R.yaml"
    script:
        "../scripts/DESeq2_transcripts.R"

rule add_gene_name:
    input:
        log = "logs/DESeq2/genes.log",
    output:
        tok = "logs/DESeq2/add_gene_name.tok"
    params:
        genes_dir = "results/DESeq2/genes_test/",
        corresponding = "data/gene_id-gene_name.txt",
        tpms = "results/kallisto/tpm.tsv"
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/add_gene_name.py"
