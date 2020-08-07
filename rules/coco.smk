rule create_coco_annotation:
    """ Create the CoCo corrected annotation """
    input:
        gtf = config["path"]['annotation']
    output:
        correct_gtf = config["path"]['correct_annotation']
    conda:
        "../envs/coco.yaml"
    shell:
        "coco correct_annotation {input.gtf}"


rule coco_correct_count:
    """ Generates the tpm counts using CoCo correct counts """
    input:
        annotation = config["path"]['correct_annotation'],
        bam_file = "results/STAR/{id}/Aligned.sortedByCoord.out.bam"
    output:
        out_file = "results/coco/{id}.tsv"
    log:
        "logs/coco/{id}.log"
    threads:
        16
    conda:
        "../envs/coco.yaml"
    shell:
        "coco cc -s 1 "
        "-t {threads} "
        "-p {input.annotation} "
        "{input.bam_file} "
        "{output.out_file}"


rule coco_merge:
    """ Merge the tpm of all conditions into one file """
    input:
        tpm_files = expand("results/coco/{id}.tsv", id=simple_id)
    output:
        counts = 'results/coco/merged/counts.tsv',
        merged = "results/coco/merged/tpm.tsv"
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/combine_coco_quantification.py"


rule coco_correct_bedgraph:
    """ Create a bedgraph from the bam files """
    input:
        bam_file = "results/STAR/{id}/Aligned.sortedByCoord.out.bam",
        chrLength = config['path']['chrNameLength']
    output:
        bedgraph = "results/coco/bedgraphs/{id}.bedgraph"
    conda:
        "../envs/coco.yaml"
    threads:
        32
    shell:
        "coco cb -u -t {threads} -c 2500000 {input.bam_file} "
        "{output.bedgraph} {input.chrLength}"


rule keep_primary_chr:
    """ Remove the first line of each bedgraph and filter out the non-standard
        chromosomes - In preparation for the bedGraphToBigWig """
    input:
        bedgraphs = expand("results/coco/bedgraphs/{id}.bedgraph", id=simple_id),
        chrLength = config['path']['chrNameLength']
    output:
        clean_bg = expand("results/coco/bedgraphs/clean_{id}.bedgraph", id=simple_id),
        new_chrLength = "data/chrNameLength_modif.txt"
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/clean_bg.py"


rule sort_bg:
    """ Sort the new bedgraphs and change the chrM to chrMT to work with bedGraphToBigWig """
    input:
        bg = "results/coco/bedgraphs/clean_{id}.bedgraph"
    output:
        clean_bg = "results/coco/bedgraphs/sorted_clean_{id}.bedgraph"
    shell:
        "sort -k1,1 -k2,2n {input.bg} | sed 's/chrM/chrMT/g' > {output.clean_bg}"


rule convert_bw:
    """ Convert the bedgraphs to bigwigs for IGV vizualisation """
    input:
        clean_bg = "results/coco/bedgraphs/sorted_clean_{id}.bedgraph",
        new_chrLength = "data/chrNameLength_modif.txt"
    output:
        bw = "results/coco/bigwig/{id}.bw"
    conda:
        "../envs/bedgraphtobigwig.yaml"
    shell:
        "bedGraphToBigWig {input.clean_bg} {input.new_chrLength} {output.bw}"
