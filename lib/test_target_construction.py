import os
import pathlib

import pytest

from lib import target_construction as tc


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
        expected = "results/fastp/PROJ1/SAM2_{}_{}_fastp.fastq".format(
            wildcards_with_lane.lane, readname
        )
    else:
        expected = "results/fastqs/PROJ1/fn4_{}.fq.gz".format(readname)
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
    expected = "fn4_{}.fq.gz".format(readname)
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
    expected = ["results/bwa-mem2/PROJ1/SAM2_{}.{}".format(x, suffix) for x in ["L001", "L002"]]
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
