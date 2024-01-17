import pandas as pd
import pytest
from snakemake.io import Namedlist


@pytest.fixture
def common_tmpdir(tmp_path):
    """
    Construct a temporary directory into which test files can be emitted
    """
    return tmp_path


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
def wildcards_local_reference():
    """
    snakemake wildcards for testing a request for a reference file with a local source
    """
    return Namedlist(fromdict={"reference_file": "verifybamid2/grch100/ref.db-UD"})


@pytest.fixture
def wildcards_s3_reference():
    """
    snakemake wildcards for testing a request for a reference file on s3
    """
    return Namedlist(fromdict={"reference_file": "deepvariant/grch100/ref.skip-regions"})


@pytest.fixture
def wildcards_url_reference():
    """
    snakemake wildcards for testing a request for a reference file specified by a URL
    """
    return Namedlist(fromdict={"reference_file": "deepvariant/grch100/ref.calling-ranges"})


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
            "snv-caller": "deepvariant",
            "sv-endpoints": {
                "strict": {
                    "sv-callers": ["manta", "lumpy", "delly"],
                    "sv-ensemble": {"min-count": 2},
                    "sv-remove-breakends": True,
                },
                "lenient": {"sv-callers": ["manta"]},
            },
            "outcome": "fastqc",
            "symlink-fastqs": True,
        },
        "references": {
            "grch100": {"fasta": "my.fa", "reportable-regions": "my.reportable-regions"}
        },
        "dnascope": {"grch100": {"model": "my.model", "dbsnp-vcf-gz": "my.vcf.gz"}},
        "verifybamid2": {
            "grch100": {"db-V": "my.V", "db-UD": "my.UD", "db-mu": "my.mu", "db-bed": "my.bed"}
        },
        "deepvariant": {
            "grch100": {
                "skip-regions": "s3://my.skip-regions",
                "calling-ranges": "https://my.calling-ranges",
            },
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
            "r1": ["fn{}_R1.fastq.gz".format(x + 1) for x in range(6)],
            "r2": ["fn{}_R2.fastq.gz".format(x + 1) for x in range(6)],
        }
    )
