import warnings

import pytest

from lib import string_manipulations as sm


@pytest.mark.parametrize(
    "input, expected",
    [
        ("grch38", "GRCh38"),
        ("GRCh38", "GRCh38"),
        ("grch37", "GRCh37"),
        ("GRCh37", "GRCh37"),
        ("hg38", "GRCh38"),
        ("HG38", "GRCh38"),
        ("hg19", "GRCh37"),
        ("HG19", "GRCh37"),
    ],
)
def test_format_reference_build(input, expected):
    """
    Test that format_reference_build can ingest a variety
    of reference genome code specifications.
    """
    output = sm.format_reference_build(input)
    assert output == expected


@pytest.mark.parametrize(
    "input, expected",
    [("hg18", "hg18"), ("grch380", "grch380")],
)
def test_format_reference_build_with_warning(input, expected):
    """
    Test that format_reference_build recapitulates its input
    when it doesn't recognize the specified build code.
    """
    with pytest.warns(UserWarning):
        output = sm.format_reference_build(input)
    assert output == expected
