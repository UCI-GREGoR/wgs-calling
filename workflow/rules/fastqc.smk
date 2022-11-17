rule run_fastqc:
    """
    Run fastQC on a paired end library
    """
    input:
        r1="results/fastqs/{prefix}_R1_{suffix}.fastq.gz",
        r2="results/fastqs/{prefix}_R2_{suffix}.fastq.gz",
    output:
        r1_zip="results/fastqc/{prefix}_R1_{suffix}_fastqc.zip",
        r2_zip="results/fastqc/{prefix}_R2_{suffix}_fastqc.zip",
        r1_html="results/fastqc/{prefix}_R1_{suffix}_fastqc.html",
        r2_html="results/fastqc/{prefix}_R2_{suffix}_fastqc.html",
    params:
        outdir="results/fastqc",
    conda:
        "../envs/fastqc.yaml"
    threads: 1
    resources:
        h_vmem="1000",
        qname="small",
    shell:
        "mkdir -p {params.outdir} && fastqc --threads {threads} {input.r1} {input.r2} --outdir {params.outdir}"
