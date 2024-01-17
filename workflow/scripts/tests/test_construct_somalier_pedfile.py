import os
import pathlib
import re
import runpy

import pandas as pd
import pytest
import snakemake.script as sms
import yaml
from snakemake.io import Namedlist


@pytest.fixture
def common_tmpdir(tmp_path):
    """
    Define a common temporary directory
    for simulating output filesystem
    """
    return tmp_path


@pytest.fixture
def ped_filename(common_tmpdir):
    """
    Name of output plink ped file for conversion function
    """
    return common_tmpdir / "construct_somalier_pedfile/output.ped"


@pytest.fixture
def problems_filename(common_tmpdir):
    """
    Name of output file reporting problematic entries during pedfile creation
    """
    return common_tmpdir / "construct_somalier_pedfile/problems.tsv"


@pytest.fixture
def linker_filename(common_tmpdir):
    """
    Name of pipeline internal linker file.
    """
    pathlib.Path(common_tmpdir / "construct_somalier_pedfile/").mkdir(parents=True, exist_ok=True)
    df = pd.DataFrame(
        {
            "subject": ["EXT{}".format(i + 1) for i in range(10)],
            "jira": ["NA" for i in range(10)],
            "project": ["proj" for i in range(10)],
            "index": ["SAM{}".format(i + 1) for i in range(10)],
            "analyte": ["analyte{}".format(i + 1) for i in range(10)],
            "sex": ["Female" if i % 2 == 0 else "Male" for i in range(10)],
        }
    )
    fn = common_tmpdir / "construct_somalier_pedfile/linker.tsv"

    df.to_csv(fn, sep="\t")
    return fn


@pytest.fixture
def standard_input(linker_filename):
    """
    Create snakemake input block containing requested linker filename
    """
    res = Namedlist(fromdict={"0": linker_filename})
    return res


@pytest.fixture
def standard_output(ped_filename, problems_filename):
    """
    Create snakemake output block containing requested ped filename
    """
    res = Namedlist(fromdict={"ped": ped_filename, "problems": problems_filename})
    return res


@pytest.fixture
def sample_list():
    """
    Create a dummy list of sample IDs
    """
    return ["SAM{}".format(i + 1) for i in range(10)]


@pytest.fixture
def standard_params(sample_list):
    """
    Create snakemake params block containing a list of sample IDs
    """
    res = Namedlist(
        fromdict={"projectid": "proj", "subjectids": sample_list, "last_sample_sex": "unknown"}
    )
    return res


@pytest.fixture
def standard_snakemake_object(standard_input, standard_output, standard_params):
    """
    Create basic snakemake object for simple runthrough test
    """
    res = sms.Snakemake(
        standard_input,
        standard_output,
        standard_params,
        Namedlist(),
        1,
        Namedlist(),
        Namedlist(),
        {},
        "",
        [],
    )
    ## override depickling
    res.input = standard_input
    res.output = standard_output
    res.params = standard_params
    return res


def test_standard_request(standard_snakemake_object):
    """
    Run simple end-to-end test with expected successful output.
    """
    pathlib.Path(os.path.dirname(standard_snakemake_object.output[0])).mkdir(
        parents=True, exist_ok=True
    )
    runpy.run_path(
        "workflow/scripts/construct_somalier_pedfile.py",
        init_globals={"snakemake": standard_snakemake_object},
    )
    assert pathlib.Path(standard_snakemake_object.output[0]).is_file()
    subjectids = standard_snakemake_object.params["subjectids"]
    expected = pd.DataFrame(
        {
            "FID": [re.sub("SAM", "EXT", i) for i in subjectids],
            "Sample": subjectids,
            "Pat": [0 for i in range(10)],
            "Mat": [0 for i in range(10)],
            "Sex": [2 if (int(re.sub("[a-zA-Z]+", "", i)) - 1) % 2 == 0 else 1 for i in subjectids],
            "Pheno": [-9 for i in range(10)],
        }
    ).set_index("FID", drop=False)
    observed = pd.read_csv(
        standard_snakemake_object.output[0],
        sep="\t",
        header=None,
        names=["FID", "Sample", "Pat", "Mat", "Sex", "Pheno"],
    ).set_index("FID", drop=False)
    assert expected.equals(observed)
