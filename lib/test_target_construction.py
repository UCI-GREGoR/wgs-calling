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
