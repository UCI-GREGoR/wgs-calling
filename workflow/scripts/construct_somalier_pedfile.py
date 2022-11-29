import pandas as pd


def run_construct_somalier_pedfile(sampleids: list, outfn: str) -> None:
    """
    Create a somalier format pedfile for linking sample id to sex annotation.

    Currently, this populates the sex entry for all subjects with placeholder 0 for unknown.
    """
    ## Due to fastqs from multiple lanes, the input sample list may contain duplicates,
    ## which is not acceptable in this context.
    sampleids_unique = list(set(sampleids))
    sampleids_unique.sort()
    x = pd.DataFrame(
        data={
            "FID": sampleids_unique,
            "Sample": sampleids_unique,
            "Pat": ["0" for x in sampleids_unique],
            "Mat": ["0" for x in sampleids_unique],
            "Sex": ["0" for x in sampleids_unique],
            "Pheno": ["-9" for x in sampleids_unique],
        }
    )
    x.to_csv(outfn, sep="\t", index=False, header=False)


run_construct_somalier_pedfile(snakemake.params["subjectids"], snakemake.output[0])  # noqa: F821
