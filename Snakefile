import os

configfile: "config.json"

original_name = list(config['datasets'].values())
simple_id = list(config['datasets'].keys())

rule all:
    input:
        counts = 'results/coco/merged/counts.tsv',
        # tok = "logs/DESeq2/add_gene_name.tok"


rule download_genome:
    """ Downloads the genome from Ensembl FTP servers """
    output:
        genome = config['path']['genome']
    params:
        link = config['download']['genome']
    shell:
        "wget --quiet -O {output.genome}.gz {params.link} && "
        "gzip -d {output.genome}.gz "


rule rename_files:
    input:
        fastq = expand("data/reads/{original_name}_R{pair}_001.fastq",
                       original_name=original_name, pair=[1, 2])
    output:
        new_name = expand("data/reads/{id}_{pair}.fastq",
                          id=simple_id, pair=[1, 2])
    run:
        for id, original in config['datasets'].items():
            for num in [1, 2]:
                old = "data/reads/{}_R{}_001.fastq".format(original, num)
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
        fasta = config["path"]["genome"],
        gtf = config["path"]['annotation']
    output:
        chrNameLength = config['path']['chrNameLength']
    params:
        dir = config['path']['star_index']
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
        bam = "results/STAR/{id}/Aligned.sortedByCoord.out.bam",
        bam_logs = "results/STAR/{id}/Log.final.out"
    params:
        index = config['path']['star_index'],
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
# include: "rules/controls.smk"

# inclue Kraken2
# include: "rules/kraken.smk"

# include DESeq
# include: "rules/DESeq2.smk"
