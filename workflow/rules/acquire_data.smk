rule copy_fastqs:
    """
    Get a local copy of fastq files before analysis
    """
    input:
        r1=lambda wildcards: tc.map_fastqs_to_manifest(wildcards, manifest, "R1"),
        r2=lambda wildcards: tc.map_fastqs_to_manifest(wildcards, manifest, "R2"),
    output:
        r1="results/fastqs/{projectid}/{prefix}_R1_{suffix}.fastq.gz",
        r2="results/fastqs/{projectid}/{prefix}_R2_{suffix}.fastq.gz",
    benchmark:
        "results/performance_benchmarks/copy_fastqs/{projectid}/{prefix}_{suffix}.fastq.tsv"
    threads: 1
    resources:
        h_vmem="500",
        qname="small",
    shell:
        "cp {input.r1} {output.r1} && cp {input.r2} {output.r2}"
