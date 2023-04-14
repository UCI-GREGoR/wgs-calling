import os
import pathlib

import pandas as pd
import pytest
from snakemake.checkpoints import Checkpoint, Checkpoints
from snakemake.io import IOFile, Namedlist, Wildcards, expand
from snakemake.rules import Rule
from snakemake.workflow import Workflow

from lib import export_data as ed


@pytest.fixture
def common_tmpdir(tmp_path):
    """
    Construct a temporary directory into which test files can be emitted
    """
    return tmp_path


@pytest.fixture
def manifest_data() -> pd.DataFrame:
    """
    Construct example manifest data in a pandas DataFrame
    """
    df = pd.DataFrame(
        {
            "projectid": ["RU00001", "RU00002", "RU00002", "RU00004"],
            "sampleid": ["internal_id2", "internal_id3", "internal_id5", "internal_id9"],
        }
    )
    return df


@pytest.fixture
def linker_data() -> pd.DataFrame:
    """
    Construct example linker data in a pandas DataFrame
    """
    df = pd.DataFrame(
        {
            "pmgrc": ["study_id1", "study_id2", "study_id3", "study_id1"],
            "jira": ["RT-0001", "RT-0001", "RT-0002", "RT-0003"],
            "ru": ["RU00001", "RU00001", "RU00002", "RU00003"],
            "sq": ["internal_id1", "internal_id2", "internal_id3", "internal_id4"],
            "ls": ["LS001", "LS002", "LS003", "LS004"],
            "sex": ["Female", "Male", "Unknown", "Female"],
            "output": ["external_id1", "external_id2", "external_id3", "external_id4"],
        }
    )
    return df


@pytest.fixture
def linker_file(common_tmpdir, linker_data: pd.DataFrame) -> pathlib.Path:
    """
    Generate a temporary file containing exemplar linker
    data, and return the temporary filename
    """
    fn = common_tmpdir / "linker_file.tsv"
    linker_data.to_csv(fn, sep="\t", index=False)
    return fn


@pytest.fixture
def linker_workflow() -> Workflow:
    """
    A snakemake rule needs to be constructed as part of an existing Workflow
    object, which accepts minimally the string filename of its snakefile.
    At least for this limited use case, the existence of the snakefile
    appears to be immaterial.
    """
    return Workflow("fake_snakefile.smk")


@pytest.fixture
def linker_rule(linker_file, linker_workflow) -> Rule:
    """
    Construct the linker file generation rule that will be flagged
    as a checkpoint to snakemake later on.

    Due to the slightly deranged recursive logic of snakemake's data structures,
    the output file's IOFile object seemingly cannot be created as a separate fixture.
    """
    my_rule = Rule("generate_linker", linker_workflow)
    my_rule._output = Namedlist(fromdict={"tsv": IOFile(str(linker_file), my_rule)})
    return my_rule


@pytest.fixture
def linker_checkpoints(linker_rule) -> Checkpoints:
    """
    Generate a snakemake checkpoints object containing sufficient
    information to extract out the result of the linker checkpoint.

    This Checkpoints object has only minimal functionality, but it should
    be able to correctly respond to the query `self.linker_rule.get().output[0]`
    """
    my_checkpoints = Checkpoints()
    my_checkpoints.register(linker_rule)
    ## override snakemake's internal safeties for determining whether
    ## the checkpoint has already been evaluated
    my_checkpoints.future_output = []
    return my_checkpoints


@pytest.fixture
def export_wildcards() -> Wildcards:
    """
    Create a dummy Wildcards object for the export data test.
    Note that the Wildcards class is just an alias for the
    standard snakemake Namedlist
    """
    return Namedlist(fromdict={"projectid": "RU00001"})


@pytest.fixture
def nonexport_wildcards() -> Wildcards:
    """
    Create a dummy Wildcards object for the nonexport data test.
    Note that the Wildcards class is just an alias for the
    standard snakemake Namedlist
    """
    return Namedlist(fromdict={"projectid": "RU00002"})


@pytest.mark.parametrize(
    "suffix,expected",
    [
        (
            "vcf.gz",
            [
                "results/export/RU00001/external_id2.vcf.gz",
            ],
        ),
        (
            "cram.md5",
            [
                "results/export/RU00001/external_id2.cram.md5",
            ],
        ),
    ],
)
def test_construct_export_files(
    export_wildcards, manifest_data, linker_checkpoints, suffix, expected
):
    """
    Test primary functionality of construct_export_files, creating
    a list of targets respecting the alignment of manifest and linker
    """
    observed = ed.construct_export_files(
        export_wildcards, manifest_data, linker_checkpoints, suffix
    )
    assert expected == observed


@pytest.mark.parametrize(
    "suffix,expected",
    [
        (
            "vcf.gz",
            [
                "results/nonexport/RU00002/internal_id5.vcf.gz",
            ],
        ),
        (
            "cram.md5",
            [
                "results/nonexport/RU00002/internal_id5.cram.md5",
            ],
        ),
    ],
)
def test_construct_nonexport_files(
    nonexport_wildcards, manifest_data, linker_checkpoints, suffix, expected
):
    """
    Test primary functionality of construct_nonexport_files, creating
    a list of targets respecting the alignment of manifest and linker
    """
    observed = ed.construct_nonexport_files(
        nonexport_wildcards, manifest_data, linker_checkpoints, suffix
    )
    assert expected == observed
