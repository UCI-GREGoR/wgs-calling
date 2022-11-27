import pandas as pd
import pytest
from snakemake.io import Namedlist


@pytest.fixture
def wildcards_with_lane():
    """
    snakemake wildcards for testing steps that require specific behaviors per-lane
    """
    return Namedlist(fromdict={"projectid": "PROJ1", "sampleid": "SAM2", "lane": "L002"})


@pytest.fixture
def wildcards_without_lane():
    """
    snakemake wildcards for testing steps that occur after lane merge
    """
    return Namedlist(fromdict={"projectid": "PROJ1", "sampleid": "SAM2"})


@pytest.fixture
def standard_config():
    """
    configuration settings imitating those of an actual workflow run
    """
    return {
        "manifest": "config/manifest.tsv",
        "multiqc-config": "config/multiqc-config.yaml",
        "genome-build": "grch100",
        "behaviors": {
            "aligner": "bwa-mem2",
            "snv-caller": "octopus",
            "sv-caller": ["manta", "tiddit"],
            "outcome": "fastqc",
            "symlink-fastqs": True,
            "trim-adapters-before-alignment": "legacy",
        },
        "references": {
            "grch100": {"fasta": "my.fa", "reportable-regions": "my.reportable-regions"}
        },
        "dnascope": {"grch100": {"model": "my.model", "dbsnp-vcf-gz": "my.vcf.gz"}},
        "verifybamid2": {
            "grch100": {"db-V": "my.V", "db-UD": "my.UD", "db-mu": "my.mu", "db-bed": "my.bed"}
        },
        "octopus": {
            "error-model": "my.error-model",
            "forest-model": "my.forest-model",
            "grch100": {"skip-regions": "my.skip-regions", "calling-ranges": "my.calling-ranges"},
        },
    }


@pytest.fixture
def standard_manifest():
    """
    Manifest table for an exemplar workflow run
    """
    return pd.DataFrame(
        {
            "projectid": ["PROJ1", "PROJ1", "PROJ1", "PROJ1", "PROJ1", "PROJ1"],
            "sampleid": ["SAM1", "SAM1", "SAM2", "SAM2", "SAM3", "SAM3"],
            "lane": ["L001", "L002", "L001", "L002", "L001", "L002"],
            "r1": ["fn{}_R1.fq.gz".format(x + 1) for x in range(6)],
            "r2": ["fn{}_R2.fq.gz".format(x + 1) for x in range(6)],
        }
    )
