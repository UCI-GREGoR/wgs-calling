## This file contains logic for merging fastqs before assessment
## with fastqc for post-alignment quality control. This action is
## required due to idiosyncrasies in how upstream users are expecting
## their output quality control data to be formatted. This action
## obfuscates all meaningful batch information by lane, and so
## is only being performed for the version of read QC that is
## present in the post-alignment multiqc report.


rule combine_input_fastqs_by_lane:
    """
    Merge input reads per-sample across lanes

    If these were reliably bgzipped, we could just cat them; but it seems like this isn't guaranteed.
    """
    input:
        lambda wildcards: tc.get_fastqs_by_lane(wildcards, checkpoints, manifest, "001.fastq.gz"),
    output:
        temp("results/fastqs_combined/pretrimming/{projectid}/{sampleid}_{readgroup}_001.fastq.gz"),
    conda:
        "../envs/bcftools.yaml" if not use_containers else None
    container:
        "{}/bcftools.sif".format(apptainer_images) if use_containers else None
    threads: config_resources["default"]["threads"]
    resources:
        mem_mb=config_resources["default"]["memory"],
        qname=lambda wildcards: rc.select_queue(
            config_resources["default"]["queue"], config_resources["queues"]
        ),
    shell:
        "cat {input} > {output}"


use rule combine_input_fastqs_by_lane as combine_fastp_fastqs_by_lane with:
    input:
        lambda wildcards: tc.get_fastqs_by_lane(wildcards, checkpoints, manifest, "fastq.fastq.gz"),
    output:
        temp("results/fastqs_combined/posttrimming/{projectid}/{sampleid}_{readgroup}_001.fastq.gz"),
