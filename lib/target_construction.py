import os

import pandas as pd
from snakemake.io import AnnotatedString, Namedlist
from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider
from snakemake.remote.S3 import RemoteProvider as S3RemoteProvider

S3 = S3RemoteProvider()
HTTP = HTTPRemoteProvider()


def map_fastq_from_project_and_sample(
    wildcards: Namedlist, config: dict, manifest: pd.DataFrame, rp: str
) -> str:
    """
    Get a particular fastq based on projectid and sampleid
    """
    ## New: allow rewiring of the DAG to provide adapter-trimmed fastqs directly to aligner
    if config["behaviors"]["trim-adapters-before-alignment"] is True:
        return "results/fastp/{}/{}_{}_{}_fastp.fastq.gz".format(
            wildcards.projectid, wildcards.sampleid, wildcards.lane, rp
        )

    query = 'projectid == "{}" and sampleid == "{}" and lane == "{}"'.format(
        wildcards.projectid, wildcards.sampleid, wildcards.lane
    )
    result = manifest.query(query)
    assert len(result) == 1
    result = result[rp.lower()]
    return "results/fastqs/{}/{}".format(wildcards.projectid, os.path.basename(result.to_list()[0]))


def map_fastqs_to_manifest(wildcards: Namedlist, manifest: pd.DataFrame, readtag: str) -> str:
    """
    Query input manifest to find path of an input fastq
    """
    query = 'projectid == "{}" and sampleid == "{}" and lane == "{}"'.format(
        wildcards.projectid, wildcards.sampleid, wildcards.lane
    )
    result = manifest.query(query)
    result = result[readtag.lower()].to_list()
    assert len(result) == 1
    return result[0]


def get_bams_by_lane(
    wildcards: Namedlist, config: dict, manifest: pd.DataFrame, suffix: str
) -> list:
    """
    For a project and sample, get all the expected bams for the subject based on manifest lanes
    """
    query = 'projectid == "{}" and sampleid == "{}"'.format(wildcards.projectid, wildcards.sampleid)
    result = manifest.query(query)
    available_lanes = result["lane"]
    result = [
        "results/aligned/{}/{}_{}.{}".format(wildcards.projectid, wildcards.sampleid, x, suffix)
        for x in available_lanes
    ]
    return result


def construct_contamination_targets(wildcards: Namedlist, manifest: pd.DataFrame) -> list:
    """
    From basic input manifest entries, construct output targets for
    a run of verifybamid2 (for contamination)
    """
    result = [
        "results/contamination/{}/{}.vb2.selfSM".format(wildcards.projectid, x)
        for x in manifest.loc[manifest["projectid"] == wildcards.projectid, "sampleid"]
    ]
    return list(set(result))


def construct_mosdepth_targets(wildcards: Namedlist, manifest: pd.DataFrame) -> list:
    """
    From basic input manifest entries, construct output targets for
    a run of mosdepth
    """
    result = [
        "results/mosdepth/{}/{}.mosdepth.global.dist.txt".format(wildcards.projectid, x)
        for x in manifest.loc[manifest["projectid"] == wildcards.projectid, "sampleid"]
    ]
    return list(set(result))


def construct_alignstats_targets(wildcards: Namedlist, manifest: pd.DataFrame) -> list:
    """
    From basic input manifest entries, construct output targets for
    a run of alignstats
    """
    result = [
        "results/alignstats/{}/{}.alignstats.json".format(wildcards.projectid, x)
        for x in manifest.loc[manifest["projectid"] == wildcards.projectid, "sampleid"]
    ]
    return list(set(result))


def construct_combined_alignstats_targets(wildcards: Namedlist) -> list:
    """
    From basic input manifest entries, construct output targets for
    combined alignstats output
    """
    ## this is somewhat idiosyncratic, but similar functions all return lists, and snakemake
    ## doesn't care, but having this return a list as well is helpful for consistency's sake.
    result = ["results/alignstats/{}/alignstats_summary_mqc.tsv".format(wildcards.projectid)]
    return result


def construct_picard_qc_targets(wildcards: Namedlist, manifest: pd.DataFrame) -> list:
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
    return list(set(result))


def construct_snv_targets(config: dict, manifest: pd.DataFrame) -> list:
    """
    From basic input configuration and manifest entries,
    construct output targets for selected caller, after
    scattergather is complete
    """
    result = [
        "results/{}/{}/{}.sorted.vcf.gz".format(config["behaviors"]["snv-caller"], x[0], x[1])
        for x in zip(manifest["projectid"], manifest["sampleid"])
    ]
    return list(set(result))


def construct_sv_targets(manifest: pd.DataFrame) -> list:
    """
    From basic input manifest entries, construct output targets for
    arbitrary SV calling
    """
    result = [
        "results/final/{}/{}.sv.vcf.gz".format(x[0], x[1])
        for x in zip(manifest["projectid"], manifest["sampleid"])
    ]
    return list(set(result))


def construct_somalier_extract_targets(wildcards: Namedlist, manifest: pd.DataFrame) -> list:
    """
    From basic input manifest entries, construct output targets for
    a run of somalier extract
    """
    result = [
        "results/somalier/{}/extract/{}.somalier".format(wildcards.projectid, x)
        for x in manifest.loc[manifest["projectid"] == wildcards.projectid, "sampleid"]
    ]
    return list(set(result))


def construct_somalier_relate_targets(wildcards: Namedlist) -> list:
    """
    From basic input manifest entries, construct output targets for
    a run of somalier relate
    """
    result = ["results/somalier/{}/relate/somalier.html".format(wildcards.projectid)]
    return list(set(result))


def construct_fastqc_targets(wildcards: Namedlist, manifest: pd.DataFrame) -> list:
    """
    From basic input manifest entries, construct output targets for
    a run of fastQC
    """
    results_prefix = "results/fastqc"
    results_r1 = [
        "{}/{}/{}_fastqc.zip".format(
            results_prefix, wildcards.projectid, os.path.basename(x).removesuffix(".fastq.gz")
        )
        for x in manifest.loc[manifest["projectid"] == wildcards.projectid, "r1"].to_list()
    ]
    results_r2 = [
        "{}/{}/{}_fastqc.zip".format(
            results_prefix, wildcards.projectid, os.path.basename(x).removesuffix(".fastq.gz")
        )
        for x in manifest.loc[manifest["projectid"] == wildcards.projectid, "r2"].to_list()
    ]
    results_r1.extend(results_r2)
    return list(set(results_r1))


def construct_fastqc_posttrimming_targets(wildcards: Namedlist, manifest: pd.DataFrame) -> list:
    """
    From basic input manifest entries, construct output targets for
    a run of fastQC, for fastqs after trimming
    """
    results_prefix = "results/fastqc_posttrimming"
    results_r1 = [
        "{}/{}/{}_fastp_fastqc.zip".format(
            results_prefix, wildcards.projectid, os.path.basename(x).removesuffix("_001.fastq.gz")
        )
        for x in manifest.loc[manifest["projectid"] == wildcards.projectid, "r1"].to_list()
    ]
    results_r2 = [
        "{}/{}/{}_fastp_fastqc.zip".format(
            results_prefix, wildcards.projectid, os.path.basename(x).removesuffix("_001.fastq.gz")
        )
        for x in manifest.loc[manifest["projectid"] == wildcards.projectid, "r2"].to_list()
    ]
    results_r1.extend(results_r2)
    return list(set(results_r1))


def construct_fastqc_combined_targets(wildcards: Namedlist, manifest: pd.DataFrame) -> list:
    """
    From basic input manifest entries, construct output targets for
    a run of fastQC, but on lane-combined fastqs
    """
    results_prefix = "results/fastqc_combined"
    results_r1 = [
        "{}/{}/{}_R1_fastqc.zip".format(results_prefix, wildcards.projectid, x)
        for x in manifest.loc[manifest["projectid"] == wildcards.projectid, "sampleid"].to_list()
    ]
    results_r2 = [
        "{}/{}/{}_R2_fastqc.zip".format(results_prefix, wildcards.projectid, x)
        for x in manifest.loc[manifest["projectid"] == wildcards.projectid, "sampleid"].to_list()
    ]
    results_r1.extend(results_r2)
    return list(set(results_r1))


def construct_fastqc_posttrimming_combined_targets(
    wildcards: Namedlist, manifest: pd.DataFrame
) -> list:
    """
    From basic input manifest entries, construct output targets for
    a run of fastQC, for fastqs after trimming, but on lane-combined fastqs
    """
    results_prefix = "results/fastqc_posttrimming_combined"
    results_r1 = [
        "{}/{}/{}_R1_fastp_fastqc.zip".format(results_prefix, wildcards.projectid, x)
        for x in manifest.loc[manifest["projectid"] == wildcards.projectid, "sampleid"].to_list()
    ]
    results_r2 = [
        "{}/{}/{}_R2_fastp_fastqc.zip".format(results_prefix, wildcards.projectid, x)
        for x in manifest.loc[manifest["projectid"] == wildcards.projectid, "sampleid"].to_list()
    ]
    results_r1.extend(results_r2)
    return list(set(results_r1))


def construct_fastp_targets(wildcards: Namedlist, manifest: pd.DataFrame) -> list:
    """
    From basic input manifest entries, construct output targets for
    a run of fastp
    """
    results_prefix = "results/fastp"
    results_r1 = [
        "{}/{}/{}_fastp.html".format(
            results_prefix,
            wildcards.projectid,
            os.path.basename(x).removesuffix(".fastq.gz").split("_R1_")[0],
        )
        for x in manifest.loc[manifest["projectid"] == wildcards.projectid, "r1"].to_list()
    ]
    results_r2 = [
        "{}/{}/{}_fastp.html".format(
            results_prefix,
            wildcards.projectid,
            os.path.basename(x).removesuffix(".fastq.gz").split("_R2_")[0],
        )
        for x in manifest.loc[manifest["projectid"] == wildcards.projectid, "r2"].to_list()
    ]
    results_r1.extend(results_r2)
    return results_r1


def map_reference_file(wildcards: Namedlist, config: dict) -> str | AnnotatedString:
    """
    Use wildcard information to figure out what configured
    reference file is needed, and then wrap that file in a
    remote provider structure as required.
    """
    queries = wildcards.reference_file.split("/")
    queries[len(queries) - 1] = queries[len(queries) - 1].removeprefix("ref.").replace(".", "-")
    current_lvl = config
    for query in queries:
        current_lvl = current_lvl[query]
    ## The intention for this function was to distinguish between S3 file paths and others,
    ## and return wrapped objects related to the remote provider service when appropriate.
    ## There have been periodic issues with the remote provider interface, but it seems
    ## to be working, somewhat inefficiently but very conveniently, for the time being.
    if current_lvl.startswith("s3://"):
        return S3.remote(current_lvl)
    elif current_lvl.startswith("https://"):
        return HTTP.remote(current_lvl)
    return current_lvl


def caller_interval_file_count(config: dict) -> int:
    """
    Get simple line count of a file; this is intended to
    count a tiny text file containing <100 interval filenames.
    """
    fn = config[config["behaviors"]["snv-caller"]][config["genome-build"]]["calling-ranges"]
    x = 0
    with open(fn, "r") as f:
        x = len(f.readlines())
    return x
