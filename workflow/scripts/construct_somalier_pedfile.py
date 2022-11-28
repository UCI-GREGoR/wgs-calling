import pandas as pd


def run_construct_somalier_pedfile(sampleids: list, outfn: str) -> None:
    """
    Create a somalier format pedfile for linking sample id to sex annotation.

    Currently, this populates the sex entry for all subjects with placeholder 0 for unknown.
    """
    x = pd.DataFrame(
        data={
            "FID": sampleids,
            "Sample": sampleids,
            "Pat": ["0" for x in sampleids],
            "Mat": ["0" for x in sampleids],
            "Sex": ["0" for x in sampleids],
            "Pheno": ["-9" for x in sampleids],
        }
    )
    x.to_csv(outfn, sep="\t", index=False, header=False)


run_construct_somalier_pedfile(snakemake.params["subjectids"], snakemake.output[0])  # noqa: F821
