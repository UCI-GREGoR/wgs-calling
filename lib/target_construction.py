import os

import pandas as pd


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
    return [results_r1, results_r2]


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
    return [results_r1, results_r2]
