import gzip
import re

import pandas as pd


def determine_variant_source(info: str) -> str:
    """
    Since apparently truvari destroys the direct information about source
    caller in a vcf record, divine the source of a variant from other
    characteristics of the record.

    As this is currently built, complete mismatches will default to manta,
    which may be undesirable if new tools are added without updating this logic.
    """
    if re.search("DELLY", info) is not None:
        return "delly"
    elif re.search("EVDNC", info) is not None:
        return "svaba"
    elif re.search("LTE=", info) is not None:
        return "tiddit"
    elif re.search("PRPOS=", info) is not None:
        return "lumpy"
    else:
        return "manta"


def process_vcf(fn: str):
    """
    Emit summary information from a vcf
    """
    tool = []
    svtype = []
    with gzip.open(fn, "rt") as f:
        for line in f.readlines():
            if line.startswith("#"):
                continue
            tokens = line.split("\t")
            variant_source = determine_variant_source(tokens[7])
            tool.append(variant_source)
            svtype_match = re.search("SVTYPE=([^;]+)", tokens[7])
            if svtype_match is not None:
                svtype.append(svtype_match[1])
            else:
                raise ValueError(
                    "cannot extract SVTYPE from vcf record info field {}".format(tokens[7])
                )
    return tool, svtype


def truvari_add_variant_sources(vcf_fn: str, collapsed_fn: str, out_fn: str) -> None:
    """
    Generate summary information about variant source data from the collapsed data vcf.
    Add that information to the true input vcf, and then emit the output to file.
    """
    tool1, svtype1 = process_vcf(vcf_fn)
    tool2, svtype2 = process_vcf(collapsed_fn)
    df = pd.DataFrame(
        {
            "svtype": [y for x in [svtype1, svtype2] for y in x],
            "tool": [y for x in [tool1, tool2] for y in x],
        }
    )
    df.to_csv(out_fn, header=False, sep="\t", index=False)


truvari_add_variant_sources(
    snakemake.input["vcf"],  # noqa: F821
    snakemake.input["collapsed"],  # noqa: F821
    snakemake.output["tsv"],  # noqa: F821
)
