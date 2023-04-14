import os
import pathlib

import pytest
from snakemake.io import AnnotatedString, Namedlist, expand
from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider
from snakemake.remote.S3 import RemoteProvider as S3RemoteProvider

from lib import target_construction as tc

S3 = S3RemoteProvider()
HTTP = HTTPRemoteProvider()


@pytest.mark.parametrize(
    "readname, use_trimmed",
    [("R1", True), ("R1", False), ("R1", "legacy"), ("R2", True), ("R2", False), ("R2", "legacy")],
)
def test_map_fastq_from_project_and_sample(
    wildcards_with_lane, standard_config, standard_manifest, readname, use_trimmed
):
    """
    Test that map_fastq_from_project_and_sample can manage
    to link between wildcards and manifest entries.
    """
    config = standard_config
    config["behaviors"]["trim-adapters-before-alignment"] = use_trimmed
    ## with trimming, the target should be fastp output; otherwise, it's the fastq input symlink or copy
    if use_trimmed is True:
        expected = "results/fastp/PROJ1/SAM2_{}_{}_fastp.fastq.gz".format(
            wildcards_with_lane.lane, readname
        )
    else:
        expected = "results/fastqs/PROJ1/fn4_{}.fastq.gz".format(readname)
    observed = tc.map_fastq_from_project_and_sample(
        wildcards_with_lane, config, standard_manifest, readname
    )
    assert observed == expected


@pytest.mark.parametrize(
    "readname",
    ["R1", "R2"],
)
def test_map_fastqs_to_manifest(
    wildcards_with_lane,
    standard_manifest,
    readname,
):
    """
    Test that map_fastqs_to_manifest can manage
    to link between wildcards and manifest entries.
    """
    expected = "fn4_{}.fastq.gz".format(readname)
    observed = tc.map_fastqs_to_manifest(wildcards_with_lane, standard_manifest, readname)
    assert observed == expected


@pytest.mark.parametrize(
    "suffix",
    ["mysuffix.bam", "mysuffix.bam.bai"],
)
def test_get_bams_by_lane(wildcards_without_lane, standard_config, standard_manifest, suffix):
    """
    Test that get_bams_by_lane can manage
    to acquire all lane-specific files for a sample
    based on manifest fastq entries.
    """
    expected = ["results/aligned/PROJ1/SAM2_{}.{}".format(x, suffix) for x in ["L001", "L002"]]
    observed = tc.get_bams_by_lane(
        wildcards_without_lane, standard_config, standard_manifest, suffix
    )
    assert observed == expected


def test_construct_contamination_targets(wildcards_without_lane, standard_manifest):
    """
    Test that construct_contamination_targets can determine
    the output files of verifybamid2.
    """
    expected = [
        "results/contamination/PROJ1/{}.vb2.selfSM".format(x) for x in ["SAM1", "SAM2", "SAM3"]
    ]
    observed = tc.construct_contamination_targets(wildcards_without_lane, standard_manifest)
    ## this function is used for snakemake target population, so order is irrelevant
    expected.sort()
    observed.sort()
    assert observed == expected


def test_construct_mosdepth_targets(wildcards_without_lane, standard_manifest):
    """
    Test that construct_mosdepth_targets can determine
    the output files of mosdepth.
    """
    expected = [
        "results/mosdepth/PROJ1/{}.mosdepth.global.dist.txt".format(x)
        for x in ["SAM1", "SAM2", "SAM3"]
    ]
    observed = tc.construct_mosdepth_targets(wildcards_without_lane, standard_manifest)
    ## this function is used for snakemake target population, so order is irrelevant
    expected.sort()
    observed.sort()
    assert observed == expected


def test_construct_alignstats_targets(wildcards_without_lane, standard_manifest):
    """
    Test that construct_alignstats_targets can determine
    the output files of alignstats.
    """
    expected = [
        "results/alignstats/PROJ1/{}.alignstats.json".format(x) for x in ["SAM1", "SAM2", "SAM3"]
    ]
    observed = tc.construct_alignstats_targets(wildcards_without_lane, standard_manifest)
    ## this function is used for snakemake target population, so order is irrelevant
    expected.sort()
    observed.sort()
    assert observed == expected


def test_construct_combined_alignstats_targets(wildcards_without_lane):
    """
    Test that construct_combined_alignstats_targets can determine
    the (post-merge) output of alignstats for a project.
    """
    expected = ["results/alignstats/PROJ1/alignstats_summary_mqc.yaml"]
    observed = tc.construct_combined_alignstats_targets(wildcards_without_lane)
    assert observed == expected


def test_construct_picard_qc_targets(wildcards_without_lane, standard_manifest):
    """
    Test that construct_picard_qc_targets can determine
    the output files for each picard post-alignment qc utility.
    """
    expected = [
        "results/collectmultiplemetrics/PROJ1/SAM1.picard.alignment_summary_metrics.txt",
        "results/collectmultiplemetrics/PROJ1/SAM2.picard.alignment_summary_metrics.txt",
        "results/collectmultiplemetrics/PROJ1/SAM3.picard.alignment_summary_metrics.txt",
        "results/collectgcbiasmetrics/PROJ1/SAM1.picard.gc_bias_metrics.txt",
        "results/collectgcbiasmetrics/PROJ1/SAM2.picard.gc_bias_metrics.txt",
        "results/collectgcbiasmetrics/PROJ1/SAM3.picard.gc_bias_metrics.txt",
        "results/collectwgsmetrics/PROJ1/SAM1.picard.collect_wgs_metrics.txt",
        "results/collectwgsmetrics/PROJ1/SAM2.picard.collect_wgs_metrics.txt",
        "results/collectwgsmetrics/PROJ1/SAM3.picard.collect_wgs_metrics.txt",
    ]
    observed = tc.construct_picard_qc_targets(wildcards_without_lane, standard_manifest)
    ## this function is used for snakemake target population, so order is irrelevant
    expected.sort()
    observed.sort()
    assert observed == expected


def test_construct_snv_targets(standard_config, standard_manifest):
    """
    Test that construct_snv_targets can determine
    the output sorted vcf files of the SNV caller sub-dag.
    """
    expected = [
        "results/deepvariant/PROJ1/{}.sorted.vcf.gz".format(x) for x in ["SAM1", "SAM2", "SAM3"]
    ]
    expected.extend(
        ["results/deepvariant/PROJ1/{}.sorted.g.vcf.gz".format(x) for x in ["SAM1", "SAM2", "SAM3"]]
    )
    observed = tc.construct_snv_targets(standard_config, standard_manifest)
    ## this function is used for snakemake target population, so order is irrelevant
    expected.sort()
    observed.sort()
    assert observed == expected


def test_construct_sv_targets(standard_manifest):
    """
    Test that construct_sv_targets can determine
    the output sorted vcf files of the structural variant sub-dag.
    """
    expected = ["results/final/PROJ1/{}.sv.vcf.gz".format(x) for x in ["SAM1", "SAM2", "SAM3"]]
    observed = tc.construct_sv_targets(standard_manifest)
    ## this function is used for snakemake target population, so order is irrelevant
    expected.sort()
    observed.sort()
    assert observed == expected


def test_construct_somalier_extract_targets(wildcards_without_lane, standard_manifest):
    """
    Test that construct_somalier_extract_targets can determine
    the output files of the first (extract) step of somalier.
    """
    expected = [
        "results/somalier/PROJ1/extract/{}.somalier".format(x) for x in ["SAM1", "SAM2", "SAM3"]
    ]
    observed = tc.construct_somalier_extract_targets(wildcards_without_lane, standard_manifest)
    ## this function is used for snakemake target population, so order is irrelevant
    expected.sort()
    observed.sort()
    assert observed == expected


def test_construct_somalier_relate_targets(wildcards_without_lane):
    """
    Test that construct_somalier_relate_targets can determine
    the output files of the second (relate) step of somalier.
    """
    expected = ["results/somalier/PROJ1/relate/somalier.html"]
    observed = tc.construct_somalier_relate_targets(wildcards_without_lane)
    ## this function is used for snakemake target population, so order is irrelevant
    expected.sort()
    observed.sort()
    assert observed == expected


def test_construct_fastqc_targets(wildcards_without_lane, standard_manifest):
    """
    Test that construct_fastqc_targets can determine
    the output files of fastqc on input R1/R2 fastq.gz files.
    """
    expected = expand(
        "results/fastqc/PROJ1/{sname}_{lname}_{readname}_001_fastqc.zip",
        sname=["SAM1", "SAM2", "SAM3"],
        lname=["L001", "L002"],
        readname=["R1", "R2"],
    )
    observed = tc.construct_fastqc_targets(wildcards_without_lane, standard_manifest)
    ## this function is used for snakemake target population, so order is irrelevant
    expected.sort()
    observed.sort()
    assert observed == expected


def test_construct_fastqc_posttrimming_targets(wildcards_without_lane, standard_manifest):
    """
    Test that construct_fastqc_posttrimming_targets can determine
    the output files of fastqc on fastp-trimmed R1/R2 fastq.gz files.
    """
    expected = expand(
        "results/fastqc_posttrimming/PROJ1/{sname}_{lname}_{readname}_fastp_fastqc.zip",
        sname=["SAM1", "SAM2", "SAM3"],
        lname=["L001", "L002"],
        readname=["R1", "R2"],
    )
    observed = tc.construct_fastqc_posttrimming_targets(wildcards_without_lane, standard_manifest)
    ## this function is used for snakemake target population, so order is irrelevant
    expected.sort()
    observed.sort()
    assert observed == expected


def test_construct_fastqc_combined_targets(wildcards_without_lane, standard_manifest):
    """
    Test that construct_fastqc_combined_targets can determine
    the output files of fastqc on input R1/R2 fastq.gz files,
    but combined across lanes.
    """
    expected = expand(
        "results/fastqc_combined/PROJ1/{sname}_{readname}_fastqc.zip",
        sname=["SAM1", "SAM2", "SAM3"],
        readname=["R1", "R2"],
    )
    observed = tc.construct_fastqc_combined_targets(wildcards_without_lane, standard_manifest)
    ## this function is used for snakemake target population, so order is irrelevant
    expected.sort()
    observed.sort()
    assert observed == expected


def test_construct_fastqc_posttrimming_combined_targets(wildcards_without_lane, standard_manifest):
    """
    Test that construct_fastqc_posttrimming_combined_targets can determine
    the output files of fastqc on fastp-trimmed R1/R2 fastq.gz files, after
    combining across lanes.
    """
    expected = expand(
        "results/fastqc_posttrimming_combined/PROJ1/{sname}_{readname}_fastqc.zip",
        sname=["SAM1", "SAM2", "SAM3"],
        readname=["R1", "R2"],
    )
    observed = tc.construct_fastqc_posttrimming_combined_targets(
        wildcards_without_lane, standard_manifest
    )
    ## this function is used for snakemake target population, so order is irrelevant
    expected.sort()
    observed.sort()
    assert observed == expected


def test_construct_fastp_targets(wildcards_without_lane, standard_manifest):
    """
    Test that construct_fastp_targets can determine
    the output files of fastp on input R1/R2 fastq.gz files.
    """
    expected = expand(
        "results/fastp/PROJ1/{sname}_{lname}_fastp.html",
        sname=["SAM1", "SAM2", "SAM3"],
        lname=["L001", "L002"],
    )
    observed = tc.construct_fastp_targets(wildcards_without_lane, standard_manifest)
    ## this function is used for snakemake target population, so order is irrelevant
    expected.sort()
    observed.sort()
    assert observed == expected


def test_map_reference_file_localfile(wildcards_local_reference, standard_config):
    """
    Test that map_reference_file can successfully detect things that look
    like local files.
    """
    expected = "my.UD"
    observed = tc.map_reference_file(wildcards_local_reference, standard_config)
    assert observed == expected


def test_map_reference_file_s3(wildcards_s3_reference, standard_config):
    """
    Test that map_reference_file can successfully detect things that look
    like they come from an s3 bucket
    """
    expected = "s3://my.skip-regions"
    observed = tc.map_reference_file(wildcards_s3_reference, standard_config)
    assert observed == expected


def test_map_reference_file_url(wildcards_url_reference, standard_config):
    """
    Test that map_reference_file can successfully detect things that look
    like they come from a URL
    """
    expected = "https://my.calling-ranges"
    observed = tc.map_reference_file(wildcards_url_reference, standard_config)
    assert observed == expected
