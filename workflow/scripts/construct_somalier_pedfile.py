import pandas as pd


def run_construct_somalier_pedfile(sampleids: list, outfn: str) -> None:
    """
    Create a somalier format pedfile for linking sample id to sex annotation.

    Currently, this populates the sex entry for all subjects with placeholder U for known.
    """
    x = pd.DataFrame(data={"Sample": sampleids, "Sex": ["U" for x in sampleids]})
    x.to_csv(outfn, sep="\t", index=False)


run_construct_somalier_pedfile(snakemake.params["subjectids"], snakemake.output[0])  # noqa: F821
