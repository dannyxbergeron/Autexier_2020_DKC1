import os

configfile: "config.json"

original_name = list(config['datasets_test'].values())
simple_id = list(config['datasets_test'].keys())

rule all:
    input:
        merged = "results/coco/merged/tpm.tsv",
        bedgraphs = expand("results/coco/bedgraphs/{id}.bedgraph",
                        id=simple_id),
        bw =  expand("results/coco/bigwig/{id}.bw",
                        id=simple_id),
        metrics = expand("results/picard/{id}/picard_insert_size_metrics.txt",
                        id=simple_id),
        merged_picard = "results/picard/picard_sumup.csv",
        merged_star = "results/STAR/star_sumup.csv",
        samtools_idx = expand("results/samtools_idxstats/{id}_idxstat.tsv",
                                id=simple_id)
                                

# rule download_genome:
#     """ Downloads the genome from Ensembl FTP servers """
#     output:
#         genome = config['path_test']['genome']
#     params:
#         link = config['download']['genome']
#     shell:
#         "wget --quiet -O {output.genome}.gz {params.link} && "
#         "gzip -d {output.genome}.gz "


rule rename_files:
    """ Rename the files with meaning names """
    input:
        fastq = expand("data/reads/{original_name}_{pair}.fastq",
                       original_name=original_name, pair=[1, 2])
    output:
        new_name = expand("data/reads/{id}_{pair}.fastq",
                          id=simple_id, pair=[1, 2])
    run:
        for id, original in config['datasets_test'].items():
            for num in [1, 2]:
                old = "data/reads/{}_{}.fastq".format(original, num)
                new_ = "data/reads/{}_{}.fastq".format(id, num)

                print(old)
                print(new_)
                os.rename(old, new_)


rule trimming:
    """ Trims the FASTQ files using Trimmomatic """
    input:
        fq1 = "data/reads/{id}_1.fastq",
        fq2 = "data/reads/{id}_2.fastq"
    output:
        fq1 = "data/trimmed/{id}_1.fastq.gz",
        fq2 = "data/trimmed/{id}_2.fastq.gz",
        unpaired_fq1 = "data/trimmed/{id}_1.unpaired.fastq.gz",
        unpaired_fq2 = "data/trimmed/{id}_2.unpaired.fastq.gz"
    params:
        options = [
            "ILLUMINACLIP:data/Adapters-PE_NextSeq.fa:2:30:10", "LEADING:5",
            "TRAILING:5", "MINLEN:45"
        ]
    log:
        "logs/trimmomatic/{id}.log"
    threads:
        32
    conda:
        "envs/trimmomatic.yaml"
    shell:
        "trimmomatic PE "
        "-threads {threads} "
        "-phred33 "
        "{input.fq1} {input.fq2} "
        "{output.fq1} {output.unpaired_fq1}  "
        "{output.fq2} {output.unpaired_fq2} "
        "{params.options} "
        "&> {log}"


rule qc:
    """ Assess the FASTQ quality using FastQC """
    input:
        fq1 = rules.trimming.output.fq1,
        fq2 = rules.trimming.output.fq2,
        unpaired_fq1 = rules.trimming.output.unpaired_fq1,
        unpaired_fq2 = rules.trimming.output.unpaired_fq2,
    output:
        fq1_out = "data/qc/{id}_1_fastqc.html"
    params:
        out_dir = "data/qc"
    log:
        "logs/fastqc/{id}.log"
    threads:
        32
    conda:
        "envs/fastqc.yaml"
    shell:
        "fastqc "
        "--outdir {params.out_dir} "
        "--format fastq "
        "--threads {threads} "
        "{input.fq1} {input.fq2} "
        "{input.unpaired_fq1} {input.unpaired_fq2} "
        "&> {log}"


rule star_index:
    """ Generates the genome index for STAR """
    input:
        fasta = config["path_test"]["genome"],
        gtf = config["path_test"]['annotation']
    output:
        chrNameLength = config['path_test']['chrNameLength']
    params:
        dir = config['path_test']['star_index']
    log:
        "logs/STAR/index.log"
    conda:
        "envs/star.yaml"
    threads:
        8
    shell:
        "mkdir -p {params.dir} && "
        "STAR --runThreadN {threads} "
        "--runMode genomeGenerate "
        "--genomeDir {params.dir} "
        "--genomeFastaFiles {input.fasta} "
        "--sjdbGTFfile {input.gtf} "
        "--sjdbOverhang 99"
        "&> {log}"


rule star_alignReads:
    """ Generates a bam file using STAR """
    input:
        idx = rules.star_index.output,
        fq1 = rules.trimming.output.fq1,
        fq2 = rules.trimming.output.fq2
    output:
        bam = "results/STAR/{id}/Aligned.sortedByCoord.out.bam"
    params:
        index = config['path_test']['star_index'],
        output_dir = "results/STAR/{id}/"
    log:
        "logs/STAR/{id}.log"
    threads:
        32
    conda:
        "envs/star.yaml"
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

# include coco
include: "rules/coco.smk"

# include quality controls
include: "rules/controls.smk"

# include DESeq
# include: "rules/DESeq2.smk"
