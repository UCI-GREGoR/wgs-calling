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
        lambda wildcards: expand(
            "results/fastqs/{{projectid}}/{{sampleid}}_{lane}_{{readgroup}}_001.fastq.gz",
            lane=manifest.loc[manifest["sampleid"] == wildcards.sampleid, "lane"],
        ),
    output:
        temp("results/fastqs_combined/pretrimming/{projectid}/{sampleid}_{readgroup}.fastq.gz"),
    conda:
        "../envs/bcftools.yaml"
    container:
        "{}/bcftools.sif".format(apptainer_images)
    threads: 1
    resources:
        mem_mb="2000",
        qname="small",
    shell:
        "gunzip -c {input} | bgzip -c > {output}"


use rule combine_input_fastqs_by_lane as combine_fastp_fastqs_by_lane with:
    input:
        lambda wildcards: expand(
            "results/fastp/{{projectid}}/{{sampleid}}_{lane}_{{readgroup}}_fastp.fastq.gz",
            lane=manifest.loc[manifest["sampleid"] == wildcards.sampleid, "lane"],
        ),
    output:
        temp("results/fastqs_combined/posttrimming/{projectid}/{sampleid}_{readgroup}.fastq.gz"),
