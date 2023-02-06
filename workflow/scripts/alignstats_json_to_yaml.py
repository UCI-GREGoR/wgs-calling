import csv
import json
import sys

import yaml


def run_alignstats_json_to_yaml(infiles, outfile):
    """
    Once upon a time, this was a python script from
    https://github.com/invitae-internal/nextflow-clo-alignstats.
    However, upon some inspection, it looks like multiqc's custom
    content support is not fully featured when reading from tsv,
    and in fact strongly prefers multiqc-specific yaml or json.
    We can provide it that.
    """
    out_dict = {}
    out_dict["id"] = "Alignstats"
    out_dict["section_nam"] = "Alignstats"
    out_dict["description"] = "Alignstats metrics."
    out_dict["format"] = "yaml"
    out_dict["plot_type"] = "table"
    out_dict["data"] = {}
    for align_json in infiles:
        with open(align_json, "r") as in_j:
            d = json.load(in_j)
            sample_name = d["InputFileName"].split("/")[-1].split(".")[0]
            ## remove InputFileName, as it does nothing and the table is already huge
            removed_key = d.pop("InputFileName", None)  # noqa: F841
            out_dict["data"][sample_name] = d
    with open(outfile, "w") as f:
        yaml.dump(out_dict, f, sort_keys=False)


run_alignstats_json_to_yaml(snakemake.input["json"], snakemake.output["yaml"])  # noqa: F821
