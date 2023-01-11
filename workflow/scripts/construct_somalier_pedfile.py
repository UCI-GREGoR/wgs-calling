import pandas as pd


def convert_sex_representation(sample_sex: str) -> int:
    """
    Take a self-reported sex representation written out as
    an English word, and convert to a plink-style integer
    representation. 0 -> Unknown, 2 -> Female, 1 -> Male
    """
    if sample_sex.lower() == "female":
        return 2
    elif sample_sex.lower() == "male":
        return 1
    else:
        return 0


def run_construct_somalier_pedfile(linker: str, ruid: str, sampleids: list, outfn: str) -> None:
    """
    Create a somalier format pedfile for linking sample id to sex annotation.

    Currently, this populates the sex entry for all subjects with placeholder 0 for unknown.
    """
    ## Due to fastqs from multiple lanes, the input sample list may contain duplicates,
    ## which is not acceptable in this context.
    ids = pd.DataFrame(
        data={"ruid": [ruid for x in range(len(sampleids))], "sampleid": sampleids}
    ).drop_duplicates()

    ## load linker information formatted from the lab logbook
    linker_data = pd.read_csv(linker, sep="\t")

    ## iterate across manifest subjects and try to find matching values
    self_reported_sex = []
    for ruid, sampleid in zip(ids["ruid"], ids["sampleid"]):
        sample_sex = linker_data.loc[
            (linker_data["ru"] == ruid) & (linker_data["sq"] == sampleid), "sex"
        ]
        if len(sample_sex) == 1:
            self_reported_sex.append(convert_sex_representation(sample_sex[0]))
        else:
            self_reported_sex.append(0)

    x = pd.DataFrame(
        data={
            "FID": ids["sampleid"],
            "Sample": ids["sampleid"],
            "Pat": ["0" for x in ids["sampleid"]],
            "Mat": ["0" for x in ids["sampleid"]],
            "Sex": self_reported_sex,
            "Pheno": ["-9" for x in ids["sampleid"]],
        }
    )
    x.to_csv(outfn, sep="\t", index=False, header=False)


run_construct_somalier_pedfile(
    snakemake.input[0],  # noqa: F821
    snakemake.params["ruid"],  # noqa: F821
    snakemake.params["subjectids"],  # noqa: F821
    snakemake.output[0],  # noqa: F821
)
