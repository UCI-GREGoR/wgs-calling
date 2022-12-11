import os
import pathlib
import runpy

import pandas
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
def standard_output(ped_filename):
    """
    Create snakemake output block containing requested ped filename
    """
    res = Namedlist(fromdict={"0": ped_filename})
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
    res = Namedlist(fromdict={"subjectids": sample_list})
    return res


@pytest.fixture
def standard_snakemake_object(standard_output, standard_params):
    """
    Create basic snakemake object for simple runthrough test
    """
    res = sms.Snakemake(
        Namedlist(),
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
    res.input = Namedlist()
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
    ## due to upstream sorting issues, this requires terrible lexical sort
    subjectids = standard_snakemake_object.params["subjectids"]
    subjectids.sort()
    expected = pandas.DataFrame(
        {
            "FID": subjectids,
            "Sample": subjectids,
            "Pat": [0 for i in range(10)],
            "Mat": [0 for i in range(10)],
            "Sex": [0 for i in range(10)],
            "Pheno": [-9 for i in range(10)],
        }
    ).set_index("FID", drop=False)
    observed = pandas.read_csv(
        standard_snakemake_object.output[0],
        sep="\t",
        header=None,
        names=["FID", "Sample", "Pat", "Mat", "Sex", "Pheno"],
    ).set_index("FID", drop=False)
    assert expected.equals(observed)
