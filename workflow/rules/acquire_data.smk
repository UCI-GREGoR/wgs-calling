rule copy_fastqs:
    """
    Get a local copy of fastq files before analysis. Note that s3 hosted inputs are
    now supported, but in that case the symlink-fastqs configuration option cannot
    be respected.
    """
    input:
        lambda wildcards: tc.map_fastqs_to_manifest(wildcards, manifest, "R" + wildcards.readgroup)
        if not tc.map_fastqs_to_manifest(
            wildcards, manifest, "R" + wildcards.readgroup
        ).startswith("s3://")
        else [],
    output:
        fastq="results/fastqs/{projectid}/{sampleid}_{lane}_R{readgroup}_{suffix}.fastq.gz",
    params:
        fastq=lambda wildcards: tc.map_fastqs_to_manifest(
            wildcards, manifest, "R" + wildcards.readgroup
        ),
        symlink_target=config["behaviors"]["symlink-fastqs"],
        profile=config["behaviors"]["import-s3"]["profile-name"]
        if "import-s3" in config["behaviors"]
        else "default",
    benchmark:
        "results/performance_benchmarks/copy_fastqs/{projectid}/{sampleid}_{lane}_R{readgroup}_{suffix}.fastq.tsv"
    threads: config_resources["default"]["threads"]
    resources:
        mem_mb=config_resources["default"]["memory"],
        qname=lambda wildcards: rc.select_queue(
            config_resources["default"]["queue"], config_resources["queues"]
        ),
    shell:
        'if [[ "{params.fastq}" == "s3://"* ]] ; then '
        "aws s3 cp --profile {params.profile} {params.fastq} {output.fastq} ; "
        'elif [[ "{params.symlink_target}" == "True" ]] ; then '
        "ln -s $(readlink -m {params.fastq}) {output.fastq} ; "
        "else cp {params.fastq} {output.fastq} ; fi"


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
