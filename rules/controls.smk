rule picard:
    """ Get the picard alignement stats from the bam files """
    input:
        bam = "results/STAR/{id}/Aligned.sortedByCoord.out.bam"
    output:
        metrics = "results/picard/{id}/picard_insert_size_metrics.txt",
        histogram = "results/picard/{id}/picard_size_histogram.pdf"
    conda:
        "../envs/picard.yaml"
    shell:
        "picard CollectInsertSizeMetrics "
        "I={input.bam} "
        "O={output.metrics} "
        "H={output.histogram}"


rule picard_sumup:
    """ Get an simple overview of the picard statistics """
    input:
        metrics = expand("results/picard/{id}/picard_insert_size_metrics.txt",
                         id=simple_id)
    output:
        merged = "results/picard/picard_sumup.csv"
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/Picard_sumup.py"


rule star_sumup:
    """ Get an simple overview of the STAR statistics """
    input:
        log_files = expand("results/STAR/{id}/Log.final.out",
                           id=simple_id)
    output:
        merged = "results/STAR/star_sumup.csv"
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/STAR_sumup.py"


rule star_alignReads_spike_in:
    """ Generates a bam to make stats for the spike in """
    input:
        fq1 = rules.trimming.output.fq1,
        fq2 = rules.trimming.output.fq2
    output:
        bam_spike = "results/STAR_spike_in/{id}/Aligned.sortedByCoord.out.bam"
    params:
        index = config['path']['star_index_spike_in'],
        output_dir = "results/STAR_spike_in/{id}/"
    log:
        "logs/STAR_spike_in/{id}.log"
    threads:
        32
    conda:
        "../envs/star.yaml"
    shell:
        "STAR --runMode alignReads "
        "--genomeDir {params.index} "
        "--readFilesIn {input.fq1} {input.fq2}  "
        "--runThreadN {threads} "
        "--readFilesCommand zcat "
        "--outReadsUnmapped Fastx "
        "--outFilterType BySJout "
        "--outStd Log "
        "--outSAMunmapped None "
        "--outSAMtype BAM SortedByCoordinate "
        "--outFileNamePrefix {params.output_dir} "
        "--outFilterScoreMinOverLread 0.3 "
        "--outFilterMatchNminOverLread 0.3 "
        "--outFilterMultimapNmax 100 "
        "--winAnchorMultimapNmax 100 "
        "--alignEndsProtrude 5 ConcordantPair "
        "&> {log}"

rule idxstats:
    """ Generates stats for the spike in """
    input:
        bam_spike = "results/STAR_spike_in/{id}/Aligned.sortedByCoord.out.bam"
    output:
        samtools_idx = "results/samtools_idxstats/{id}_idxstat.tsv"
    conda:
        "../envs/coco.yaml"
    shell:
        "samtools index {input.bam_spike} && "
        "samtools idxstats {input.bam_spike} > {output.samtools_idx}"
