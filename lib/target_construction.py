import os

import pandas as pd
from snakemake.remote.S3 import RemoteProvider as S3RemoteProvider

S3 = S3RemoteProvider()


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


def construct_contamination_targets(wildcards, manifest: pd.DataFrame) -> list:
    """
    From basic input manifest entries, construct output targets for
    a run of verifybamid2 (for contamination)
    """
    result = [
        "results/contamination/{}/{}.vb2.selfSM".format(wildcards.projectid, x)
        for x in manifest.loc[manifest["projectid"] == wildcards.projectid, "sampleid"]
    ]
    return result


def construct_mosdepth_targets(wildcards, manifest: pd.DataFrame) -> list:
    """
    From basic input manifest entries, construct output targets for
    a run of mosdepth
    """
    result = [
        "results/mosdepth/{}/{}.mosdepth.global.dist.txt".format(wildcards.projectid, x)
        for x in manifest.loc[manifest["projectid"] == wildcards.projectid, "sampleid"]
    ]
    return result


def construct_alignstats_targets(wildcards, manifest: pd.DataFrame) -> list:
    """
    From basic input manifest entries, construct output targets for
    a run of alignstats
    """
    result = [
        "results/alignstats/{}/{}.bwa2a.alignstats.json".format(wildcards.projectid, x)
        for x in manifest.loc[manifest["projectid"] == wildcards.projectid, "sampleid"]
    ]
    return result


def construct_combined_alignstats_targets(wildcards) -> list:
    """
    From basic input manifest entries, construct output targets for
    combined alignstats output
    """
    result = ["results/alignstats/{}/alignstats_summary_mqc.tsv".format(wildcards.projectid)]
    return result


def construct_picard_qc_targets(wildcards, manifest: pd.DataFrame) -> list:
    """
    From basic input manifest entries, construct output targets for
    various QC passes with picard
    """
    result1 = [
        "results/collectmultiplemetrics/{}/{}.picard.alignment_summary_metrics.txt".format(
            wildcards.projectid, x
        )
        for x in manifest.loc[manifest["projectid"] == wildcards.projectid, "sampleid"]
    ]
    result2 = [
        "results/collectgcbiasmetrics/{}/{}.picard.gc_bias_metrics.txt".format(
            wildcards.projectid, x
        )
        for x in manifest.loc[manifest["projectid"] == wildcards.projectid, "sampleid"]
    ]
    result3 = [
        "results/collectwgsmetrics/{}/{}.picard.collect_wgs_metrics.txt".format(
            wildcards.projectid, x
        )
        for x in manifest.loc[manifest["projectid"] == wildcards.projectid, "sampleid"]
    ]
    result = [result1, result2, result3]
    result = [y for x in result for y in x]
    return result


def construct_somalier_extract_targets(wildcards, manifest: pd.DataFrame) -> list:
    """
    From basic input manifest entries, construct output targets for
    a run of somalier extract
    """
    result = [
        "results/somalier/{}/extract/{}.somalier".format(wildcards.projectid, x)
        for x in manifest.loc[manifest["projectid"] == wildcards.projectid, "sampleid"]
    ]
    return result


def construct_somalier_relate_targets(wildcards) -> list:
    """
    From basic input manifest entries, construct output targets for
    a run of somalier relate
    """
    result = ["results/somalier/{}/relate/somalier.html".format(wildcards.projectid)]
    return result


def construct_fastqc_targets(wildcards, manifest: pd.DataFrame) -> list:
    """
    From basic input manifest entries, construct output targets for
    a run of fastQC
    """
    results_prefix = "results/fastqc"
    results_r1 = [
        "{}/{}/{}_fastqc.zip".format(
            results_prefix, wildcards.projectid, os.path.basename(x).rstrip(".fastq.gz")
        )
        for x in manifest.loc[manifest["projectid"] == wildcards.projectid, "r1"].to_list()
    ]
    results_r2 = [
        "{}/{}/{}_fastqc.zip".format(
            results_prefix, wildcards.projectid, os.path.basename(x).rstrip(".fastq.gz")
        )
        for x in manifest.loc[manifest["projectid"] == wildcards.projectid, "r2"].to_list()
    ]
    results_r1.extend(results_r2)
    return results_r1


def construct_fastp_targets(wildcards, manifest: pd.DataFrame) -> list:
    """
    From basic input manifest entries, construct output targets for
    a run of fastp
    """
    results_prefix = "results/fastp"
    results_r1 = [
        "{}/{}/{}_fastp.html".format(
            results_prefix,
            wildcards.projectid,
            os.path.basename(x).rstrip(".fastq.gz").split("_R1_")[0],
        )
        for x in manifest.loc[manifest["projectid"] == wildcards.projectid, "r1"].to_list()
    ]
    results_r2 = [
        "{}/{}/{}_fastp.html".format(
            results_prefix,
            wildcards.projectid,
            os.path.basename(x).rstrip(".fastq.gz").split("_R2_")[0],
        )
        for x in manifest.loc[manifest["projectid"] == wildcards.projectid, "r2"].to_list()
    ]
    results_r1.extend(results_r2)
    return results_r1


def map_reference_file(wildcards, config: dict):
    """
    Use wildcard information to figure out what configured
    reference file is needed, and then wrap that file in a
    remote provider structure as required.
    """
    queries = wildcards.reference_file.split("/")
    queries[len(queries) - 1] = queries[len(queries) - 1].replace(".", "-")
    current_lvl = config
    for query in queries:
        current_lvl = current_lvl[query]
    if current_lvl.startswith("s3://"):
        return S3.remote(current_lvl)
    return current_lvl
