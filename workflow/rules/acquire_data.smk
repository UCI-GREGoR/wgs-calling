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
    threads: config_resources["default"]["threads"]
    resources:
        mem_mb=config_resources["default"]["memory"],
        qname=lambda wildcards: rc.select_queue(
            config_resources["default"]["queue"], config_resources["queues"]
        ),
    shell:
        'if [[ "{params.symlink_target}" == "True" ]] ; then '
        "ln -s $(readlink -m {input.r1}) {output.r1} && ln -s $(readlink -m {input.r2}) {output.r2} ; "
        "else cp {input.r1} {output.r1} && cp {input.r2} {output.r2} ; fi"


rule copy_bams:
    """
    When input files are bams, grab copies of them from external sources. Unlike with fastqs,
    these are temp flagged in favor of keeping their downstream converted fastqs.
    """
    output:
        bam=temp("results/imported_bams/{projectid}/{sampleid}.bam"),
    params:
        bam=lambda wildcards: tc.locate_input_bam(wildcards, manifest, False),
        profile=config["behaviors"]["import-s3"]["profile-name"]
        if "import-s3" in config["behaviors"]
        else "default",
    benchmark:
        "results/performance_benchmarks/copy_bams/{projectid}/{sampleid}.tsv"
    conda:
        "../envs/awscli.yaml" if not use_containers else None
    container:
        "{}/awscli.sif".format(apptainer_images) if use_containers else None
    threads: config_resources["default"]["threads"]
    resources:
        mem_mb=config_resources["default"]["memory"],
        qname=lambda wildcards: rc.select_queue(
            config_resources["default"]["queue"], config_resources["queues"]
        ),
    shell:
        'if [[ "{params.bam}" == "s3://"* ]] ; then '
        "aws s3 cp --profile {params.profile} {params.bam} {output.bam} ; "
        "else cp {params.bam} {output.bam} ; "
        "fi"
