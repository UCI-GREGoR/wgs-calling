rule run_fastqc_pretrimming:
    """
    Run fastQC on a paired end library, before trimming
    """
    input:
        r1="results/fastqs/{projectid}/{prefix}_R1_{suffix}.fastq.gz",
        r2="results/fastqs/{projectid}/{prefix}_R2_{suffix}.fastq.gz",
    output:
        r1_zip="results/fastqc/{projectid}/{prefix}_R1_{suffix}_fastqc.zip",
        r2_zip="results/fastqc/{projectid}/{prefix}_R2_{suffix}_fastqc.zip",
        r1_html="results/fastqc/{projectid}/{prefix}_R1_{suffix}_fastqc.html",
        r2_html="results/fastqc/{projectid}/{prefix}_R2_{suffix}_fastqc.html",
    benchmark:
        "results/performance_benchmarks/run_fastqc_pretrimming/{projectid}/{prefix}_{suffix}_fastqc.tsv"
    params:
        outdir="results/fastqc/{projectid}",
        tmpdir="temp",
    conda:
        "../envs/fastqc.yaml" if not use_containers else None
    container:
        "{}/fastqc.sif".format(apptainer_images) if use_containers else None
    threads: config_resources["fastqc"]["threads"]
    resources:
        mem_mb=config_resources["fastqc"]["memory"],
        qname=lambda wildcards: rc.select_queue(
            config_resources["fastqc"]["queue"], config_resources["queues"]
        ),
    shell:
        "mkdir -p {params.outdir} && "
        "mkdir -p {params.tmpdir} && "
        "fastqc --threads {threads} --outdir {params.outdir} --dir {params.tmpdir} {input.r1} {input.r2}"


rule run_fastqc_posttrimming:
    """
    Run fastQC on a paired end library, after trimming
    """
    input:
        r1="results/fastp/{projectid}/{sampleid}_{lane}_R1_fastp.fastq.gz",
        r2="results/fastp/{projectid}/{sampleid}_{lane}_R2_fastp.fastq.gz",
    output:
        r1_zip="results/fastqc_posttrimming/{projectid}/{sampleid}_{lane}_R1_fastp_fastqc.zip",
        r2_zip="results/fastqc_posttrimming/{projectid}/{sampleid}_{lane}_R2_fastp_fastqc.zip",
        r1_html="results/fastqc_posttrimming/{projectid}/{sampleid}_{lane}_R1_fastp_fastqc.html",
        r2_html="results/fastqc_posttrimming/{projectid}/{sampleid}_{lane}_R2_fastp_fastqc.html",
    benchmark:
        "results/performance_benchmarks/run_fastqc_posttrimming/{projectid}/{sampleid}_{lane}_fastp_fastqc.tsv"
    params:
        outdir="results/fastqc_posttrimming/{projectid}",
        tmpdir="temp",
    conda:
        "../envs/fastqc.yaml" if not use_containers else None
    container:
        "{}/fastqc.sif".format(apptainer_images) if use_containers else None
    threads: config_resources["fastqc"]["threads"]
    resources:
        mem_mb=config_resources["fastqc"]["memory"],
        qname=lambda wildcards: rc.select_queue(
            config_resources["fastqc"]["queue"], config_resources["queues"]
        ),
    shell:
        "mkdir -p {params.outdir} && "
        "mkdir -p {params.tmpdir} && "
        "fastqc --threads {threads} --outdir {params.outdir} --dir {params.tmpdir} {input.r1} {input.r2}"


use rule run_fastqc_pretrimming as run_fastqc_pretrimming_combined with:
    input:
        r1="results/fastqs_combined/pretrimming/{projectid}/{sampleid}_R1.fastq.gz",
        r2="results/fastqs_combined/pretrimming/{projectid}/{sampleid}_R2.fastq.gz",
    output:
        r1_zip="results/fastqc_combined/{projectid}/{sampleid}_R1_fastqc.zip",
        r2_zip="results/fastqc_combined/{projectid}/{sampleid}_R2_fastqc.zip",
        r1_html="results/fastqc_combined/{projectid}/{sampleid}_R1_fastqc.html",
        r2_html="results/fastqc_combined/{projectid}/{sampleid}_R2_fastqc.html",
    benchmark:
        "results/performance_benchmarks/run_fastqc_pretrimming_combined/{projectid}/{sampleid}_fastqc.tsv"
    params:
        outdir="results/fastqc_combined/{projectid}",
        tmpdir="tmp",


use rule run_fastqc_posttrimming as run_fastqc_posttrimming_combined with:
    input:
        r1="results/fastqs_combined/posttrimming/{projectid}/{sampleid}_R1.fastq.gz",
        r2="results/fastqs_combined/posttrimming/{projectid}/{sampleid}_R2.fastq.gz",
    output:
        r1_zip="results/fastqc_posttrimming_combined/{projectid}/{sampleid}_R1_fastqc.zip",
        r2_zip="results/fastqc_posttrimming_combined/{projectid}/{sampleid}_R2_fastqc.zip",
        r1_html="results/fastqc_posttrimming_combined/{projectid}/{sampleid}_R1_fastqc.html",
        r2_html="results/fastqc_posttrimming_combined/{projectid}/{sampleid}_R2_fastqc.html",
    benchmark:
        "results/performance_benchmarks/run_fastqc_posttrimming_combined/{projectid}/{sampleid}_fastqc.tsv"
    params:
        outdir="results/fastqc_posttrimming_combined/{projectid}",
        tmpdir="tmp",
