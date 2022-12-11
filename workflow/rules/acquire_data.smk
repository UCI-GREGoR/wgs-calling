rule copy_fastqs:
    """
    Get a local copy of fastq files before analysis
    """
    input:
        r1=lambda wildcards: tc.map_fastqs_to_manifest(wildcards, manifest, "R1"),
        r2=lambda wildcards: tc.map_fastqs_to_manifest(wildcards, manifest, "R2"),
    output:
        r1="results/fastqs/{projectid}/{sampleid}_{lane}_R1_{suffix}.fastq.gz",
        r2="results/fastqs/{projectid}/{sampleid}_{lane}_R2_{suffix}.fastq.gz",
    params:
        symlink_target=config["behaviors"]["symlink-fastqs"],
    benchmark:
        "results/performance_benchmarks/copy_fastqs/{projectid}/{sampleid}_{lane}_{suffix}.fastq.tsv"
    threads: 1
    resources:
        mem_mb="500",
        qname="small",
    shell:
        'if [[ "{params.symlink_target}" == "True" ]] ; then '
        "ln -s $(readlink -m {input.r1}) {output.r1} && ln -s $(readlink -m {input.r2}) {output.r2} ; "
        "else cp {input.r1} {output.r1} && cp {input.r2} {output.r2} ; fi"
