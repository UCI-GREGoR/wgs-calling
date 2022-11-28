import csv
import json
import sys
from pprint import pprint


def run_alignstats_json_to_tsv(infiles, outfile):
    """
    Port python script from https://github.com/invitae-internal/nextflow-clo-alignstats
    with the goal of eventually writing tests for it.
    """
    FIRST = True

    with open(outfile, "w") as out_f:
        w = csv.writer(out_f, delimiter="\t")
        for align_json in infiles:
            with open(align_json, "r") as in_j:
                d = json.load(in_j)
                samp = d["InputFileName"].split("/")[-1].split(".")[0]
                if FIRST:
                    # Write header
                    w.writerow(["# id: 'Alignstats'"])
                    w.writerow(["# section_name: 'Alignstats'"])
                    w.writerow(["# description: 'Alignstats metrics.'"])
                    w.writerow(["# format: 'tsv'"])
                    w.writerow(["# plot_type: 'table'"])
                    w.writerow(["sample"] + list(d.keys()))
                    FIRST = False
                w.writerow([samp] + list(d.values()))


run_alignstats_json_to_tsv(snakemake.input["json"], snakemake.output["tsv"])  # noqa: F821
