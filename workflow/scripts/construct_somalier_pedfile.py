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


def add_problem(problems: dict, sampleid: str, problem: str) -> dict:
    """
    Add annotation of problematic ID/annotation configuration to
    aggregator for later reporting
    """
    if sampleid in problems.keys():
        problems[sampleid] = problems[sampleid] + "; " + problem
    else:
        problems[sampleid] = problem
    return problems


def run_construct_somalier_pedfile(
    linker: str, ruid: str, sampleids: list, last_sample_sex: str, outfn: str, problemfn: str
) -> None:
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

    ## track detected problems with ID consistency
    problems = {}

    ## iterate across manifest subjects and try to find matching values
    self_reported_sex = []
    family_id = []
    mat_id = []
    pat_id = []
    parent_data = {}
    for subjectid, ruid, sqid in zip(linker_data["subject"], linker_data["ru"], linker_data["sq"]):
        parsed_sample_id = subjectid.split("-")
        ## new: permit non-compliant IDs to not actually break everything
        if len(parsed_sample_id) == 4:
            parent_data["{}_{}-{}".format(ruid, parsed_sample_id[2], parsed_sample_id[3])] = sqid

    for ruid, sampleid in zip(ids["ruid"], ids["sampleid"]):
        sample_sex = linker_data.loc[
            (linker_data["ru"] == ruid) & (linker_data["sq"] == sampleid), "sex"
        ]
        subject_id = linker_data.loc[
            (linker_data["ru"] == ruid) & (linker_data["sq"] == sampleid), "subject"
        ]
        if len(sample_sex) == 1:
            ## new: deal with IDs that don't actually look like "PMGRC-\d+-\d+-\d"
            if subject_id.to_list()[0].startswith("PMGRC"):
                parsed_sample_id = subject_id.to_list()[0].split("-")
            else:
                parsed_sample_id = [
                    "notPMGRC",
                    subject_id.to_list()[0],
                    subject_id.to_list()[0],
                    "0",
                ]
            invalid_family_structure = False
            if parsed_sample_id[1] == parsed_sample_id[2] and parsed_sample_id[3] != "0":
                problems = add_problem(
                    problems, sampleid, "subject with proband ID is not flagged -0"
                )
                invalid_family_structure = True
            sex_representation = convert_sex_representation(sample_sex.to_list()[0])
            invalid_sex_configuration = False
            if sex_representation != 0:
                if parsed_sample_id[3] == "1" and sex_representation != 1:
                    problems = add_problem(
                        problems,
                        sampleid,
                        "subject is proband father without self-reported male sex",
                    )
                    invalid_sex_configuration = True
                elif parsed_sample_id[3] == "2" and sex_representation != 2:
                    problems = add_problem(
                        problems,
                        sampleid,
                        "subject is proband mother without self-reported female sex",
                    )
                    invalid_sex_configuration = True
                if not invalid_sex_configuration:
                    self_reported_sex.append(sex_representation)
                else:
                    self_reported_sex.append(0)
            elif parsed_sample_id[3] == "1":
                self_reported_sex.append(convert_sex_representation("Male"))
            elif parsed_sample_id[3] == "2":
                self_reported_sex.append(convert_sex_representation("Female"))
            else:
                self_reported_sex.append(0)
            if (
                "{}_{}-1".format(ruid, parsed_sample_id[1]) in parent_data
                and not invalid_family_structure
            ):
                pat_id.append(parent_data["{}_{}-1".format(ruid, parsed_sample_id[1])])
            else:
                pat_id.append("0")
            if (
                "{}_{}-2".format(ruid, parsed_sample_id[1]) in parent_data
                and not invalid_family_structure
            ):
                mat_id.append(parent_data["{}_{}-2".format(ruid, parsed_sample_id[1])])
            else:
                mat_id.append("0")
            family_id.append(parsed_sample_id[2])
        else:
            self_reported_sex.append(0)
            mat_id.append(0)
            pat_id.append(0)
            family_id.append(sampleid)
            problems = add_problem(problems, sampleid, "self-reported sex missing from annotations")

    last_sample_sex = convert_sex_representation(last_sample_sex)
    if last_sample_sex != 0:
        self_reported_sex[-1] = last_sample_sex

    x = pd.DataFrame(
        data={
            "FID": family_id,
            "Sample": ids["sampleid"],
            "Pat": pat_id,
            "Mat": mat_id,
            "Sex": self_reported_sex,
            "Pheno": ["-9" for x in ids["sampleid"]],
        }
    )
    x.to_csv(outfn, sep="\t", index=False, header=False)
    for problem_key in problems.keys():
        problems[problem_key] = [problems[problem_key]]
    problems = pd.DataFrame(data=problems)
    problems.to_csv(problemfn, sep="\t", index=False, header=True)


run_construct_somalier_pedfile(
    snakemake.input[0],  # noqa: F821
    snakemake.params["ruid"],  # noqa: F821
    snakemake.params["subjectids"],  # noqa: F821
    snakemake.params["last_sample_sex"],  # noqa: F821
    snakemake.output["ped"],  # noqa: F821
    snakemake.output["problems"],  # noqa: F821
)
