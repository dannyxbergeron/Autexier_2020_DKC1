rule DESeq2_genes:
    """ Differential expression for the different conditions """
    input:
        counts = 'results/coco/merged/counts.tsv',
        samples = "data/design.tsv"
    output:
        results = directory("results/DESeq2/genes"),
    log:
        "logs/DESeq2/genes.log"
    conda:
        "../envs/R.yaml"
    script:
        "../scripts/DESeq2_genes.R"


rule add_gene_name:
    input:
        log = "logs/DESeq2/genes.log",
    output:
        tok = "logs/DESeq2/add_gene_name.tok"
    params:
        genes_dir = "results/DESeq2/genes/",
        corresponding = "data/gene_id-gene_name.txt",
        tpms = 'results/coco/merged/id_tpm.tsv'
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/add_gene_name.py"
