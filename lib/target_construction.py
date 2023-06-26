import os

import pandas as pd
from snakemake.checkpoints import Checkpoint, Checkpoints
from snakemake.io import AnnotatedString, Namedlist, expand


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
    query = 'projectid == "{}" and sampleid == "{}"'.format(wildcards.projectid, wildcards.sampleid)
    result = manifest.query(query)
    ## try to determine if one of a series of special data types are present
    ## "bam": data were input as pre-aligned bam that needs to be realigned
    if "bam" in manifest.columns:
        return "results/fastqs_from_bams/{}_{}_{}_001.fastq.gz".format(
            wildcards.sampleid, wildcards.lane, readtag
        )

    result_bylane = result.loc[result["lane"] == wildcards.lane,]
    if len(result_bylane) == 1:
        return result_bylane[readtag.lower()].to_list()[0]
    if len(result.loc[result["lane"] == "combined",]) == 1:
        return "results/fastqs_from_fastq/{}/{}_{}_{}_001.fastq.gz".format(
            wildcards.projectid, wildcards.sampleid, wildcards.lane, readtag
        )
    raise ValueError(
        "no valid manifest entry found for project {}, sample {}, lane {}, read {}".format(
            wildcards.projectid, wildcards.sampleid, wildcards.lane, readtag
        )
    )


def locate_input_bam(wildcards: Namedlist, manifest: pd.DataFrame):
    """
    Attempt to find user-specified bamfile for input. This is expected
    to be specified as an alternative to fastqs; fastqs are preferred.
    """
    query = 'projectid == "{}" and sampleid == "{}"'.format(wildcards.projectid, wildcards.sampleid)
    result = manifest.query(query)
    assert not result["bam"].isna().to_list()[0]
    return result["bam"].to_list()[0]


def get_bams_by_lane(
    wildcards: Namedlist,
    checkpoints: Checkpoints,
    config: dict,
    manifest: pd.DataFrame,
    suffix: str,
) -> list:
    """
    For a project and sample, get all the expected bams for the subject based on manifest lanes
    """
    query = 'projectid == "{}" and sampleid == "{}"'.format(wildcards.projectid, wildcards.sampleid)
    result = manifest.query(query)
    if "bam" in manifest.columns:
        with checkpoints.input_bam_sample_lanes.get(
            projectid=wildcards.projectid, sampleid=wildcards.sampleid
        ).output[0].open() as f:
            available_lanes = ["L" + x.rstrip().zfill(3) for x in f.readlines()]
    else:
        available_lanes = result["lane"]
        if len(available_lanes) == 1:
            ## determine if one of a series of special lane types is present.
            ## if any of them are present, will probably need a checkpoint's output
            ## to determine the expected set of lanes.
            available_lane = available_lanes.to_list()[0]
            if available_lane == "combined":
                with checkpoints.input_fastq_sample_lanes.get(
                    projectid=wildcards.projectid, sampleid=wildcards.sampleid, readgroup="R1"
                ).output[0].open() as f:
                    available_lanes = ["L" + x.rstrip().zfill(3) for x in f.readlines()]
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
    result = ["results/alignstats/{}/alignstats_summary_mqc.yaml".format(wildcards.projectid)]
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
    result.extend(
        [
            "results/{}/{}/{}.sorted.g.vcf.gz".format(config["behaviors"]["snv-caller"], x[0], x[1])
            for x in zip(manifest["projectid"], manifest["sampleid"])
        ]
    )
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


def construct_fastqc_targets(
    wildcards: Namedlist,
    manifest: pd.DataFrame,
    checkpoints: Checkpoints,
    results_prefix: str,
    results_suffix: str,
    include_lane: bool,
) -> list:
    """
    From basic input manifest entries, construct output targets for
    a run of fastQC. Depending on user configuration, mess around with output
    paths to make unique results combined by lane, using output from fastp,
    or both.
    """
    results = []
    for projectid, sampleid, lane in zip(
        manifest["projectid"], manifest["sampleid"], manifest["lane"]
    ):
        if projectid == wildcards.projectid:
            if include_lane:
                available_lanes = lane
                if "bam" in manifest.columns:
                    with checkpoints.input_bam_sample_lanes.get(
                        projectid=projectid, sampleid=sampleid
                    ).output[0].open() as f:
                        available_lanes = ["L" + x.rstrip().zfill(3) for x in f.readlines()]
                elif available_lanes == "combined":
                    with checkpoints.input_fastq_sample_lanes.get(
                        projectid=projectid, sampleid=sampleid, readgroup="R1"
                    ).output[0].open() as f:
                        available_lanes = ["L" + x.rstrip().zfill(3) for x in f.readlines()]
                results.extend(
                    expand(
                        "{resultdir}/{projectdir}/{sampleid}_{lane}_{readgroup}_{suffix}.zip",
                        resultdir=results_prefix,
                        projectdir=projectid,
                        sampleid=sampleid,
                        lane=available_lanes,
                        readgroup=["R1", "R2"],
                        suffix=results_suffix,
                    )
                )
            else:
                results.extend(
                    expand(
                        "{resultdir}/{projectdir}/{sampleid}_{readgroup}_{suffix}.zip",
                        resultdir=results_prefix,
                        projectdir=projectid,
                        sampleid=sampleid,
                        readgroup=["R1", "R2"],
                        suffix=results_suffix,
                    )
                )
    return list(set(results))


def construct_fastp_targets(
    wildcards: Namedlist, manifest: pd.DataFrame, checkpoints: Checkpoints
) -> list:
    """
    From basic input manifest entries, construct output targets for
    a run of fastp
    """
    results_prefix = "results/fastp"
    results = []
    for projectid, sampleid, lane in zip(
        manifest["projectid"], manifest["sampleid"], manifest["lane"]
    ):
        if projectid == wildcards.projectid:
            available_lanes = lane
            if "bam" in manifest.columns:
                with checkpoints.input_bam_sample_lanes.get(
                    projectid=projectid, sampleid=sampleid
                ).output[0].open() as f:
                    available_lanes = ["L" + x.rstrip().zfill(3) for x in f.readlines()]
            elif available_lanes == "combined":
                with checkpoints.input_fastq_sample_lanes.get(
                    projectid=projectid, sampleid=sampleid, readgroup="R1"
                ).output[0].open() as f:
                    available_lanes = ["L" + x.rstrip().zfill(3) for x in f.readlines()]
            results.extend(
                expand(
                    "{resultdir}/{projectdir}/{sampleid}_{lane}_fastp.html",
                    resultdir=results_prefix,
                    projectdir=projectid,
                    sampleid=sampleid,
                    lane=available_lanes,
                )
            )
    return list(set(results))


def map_reference_file(wildcards: Namedlist, config: dict):
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
    ## There have been periodic issues with the remote provider interface, and the FTP one
    ## is completely unusable with conda env creation due to timeouts, so I'm giving up
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
