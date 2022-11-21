import os

import pandas as pd


def map_fastq_from_project_and_sample(wildcards, manifest, rp) -> str:
    """
    Get a particular fastq based on projectid and sampleid
    """
    query = 'projectid == "{}" and sampleid == "{}"'.format(wildcards.projectid, wildcards.sampleid)
    result = manifest.query(query)
    assert len(result) == 1
    if rp == "R1":
        result = result["r1"]
    else:
        result = result["r2"]
    return "results/fastqs/{}/{}".format(wildcards.projectid, os.path.basename(result.to_list()[0]))


def map_fastqs_to_sampleid(wildcards) -> str:
    """
    Determine sample ID from fastq filename
    """
    return os.path.basename(wildcards.prefix).split("_L0")[0]


def map_fastqs_to_manifest(wildcards, manifest, readtag) -> str:
    """
    Query input manifest to find path of an input fastq
    """
    result = ""
    query = "{}_{}_{}.fastq.gz".format(wildcards.prefix, readtag, wildcards.suffix)
    if readtag == "R1":
        result = [
            x[1]
            for x in zip(manifest["projectid"], manifest["r1"])
            if x[0] == wildcards.projectid and os.path.basename(x[1]) == query
        ]
    else:
        result = [
            x[1]
            for x in zip(manifest["projectid"], manifest["r2"])
            if x[0] == wildcards.projectid and os.path.basename(x[1]) == query
        ]
    assert len(result) == 1
    return result[0]


def construct_bwamem2_targets(manifest: pd.DataFrame) -> list:
    """
    From basic input manifest entries, construct output targets for
    a run of bwa-mem2
    """
    result = [
        "results/bwa-mem2/{}/{}.bwa2a.bam".format(x[0], x[1])
        for x in zip(manifest["projectid"], manifest["sampleid"])
    ]
    return result


def construct_markdups_targets(manifest: pd.DataFrame) -> list:
    """
    From basic input manifest entries, construct output targets for
    a run of samtools markdups
    """
    result = [
        "results/markdups/{}/{}.mrkdup.sort.bam".format(x[0], x[1])
        for x in zip(manifest["projectid"], manifest["sampleid"])
    ]
    return result


def construct_fastqc_targets(manifest: pd.DataFrame) -> list:
    """
    From basic input manifest entries, construct output targets for
    a run of fastQC
    """
    results_prefix = "results/fastqc"
    results_r1 = [
        "{}/{}/{}_fastqc.zip".format(
            results_prefix, x[0], os.path.basename(x[1]).rstrip(".fastq.gz")
        )
        for x in zip(manifest["projectid"], manifest["r1"])
    ]
    results_r2 = [
        "{}/{}/{}_fastqc.zip".format(
            results_prefix, x[0], os.path.basename(x[1]).rstrip(".fastq.gz")
        )
        for x in zip(manifest["projectid"], manifest["r2"])
    ]
    results_r1.extend(results_r2)
    return results_r1


def construct_fastp_targets(manifest: pd.DataFrame) -> list:
    """
    From basic input manifest entries, construct output targets for
    a run of fastp
    """
    results_prefix = "results/fastp"
    results_r1 = [
        "{}/{}/{}_fastp.html".format(
            results_prefix, x[0], os.path.basename(x[1]).rstrip(".fastq.gz").split("_R1_")[0]
        )
        for x in zip(manifest["projectid"], manifest["r1"])
    ]
    results_r2 = [
        "{}/{}/{}_fastp.html".format(
            results_prefix, x[0], os.path.basename(x[1]).rstrip(".fastq.gz").split("_R2_")[0]
        )
        for x in zip(manifest["projectid"], manifest["r2"])
    ]
    results_r1.extend(results_r2)
    return results_r1
