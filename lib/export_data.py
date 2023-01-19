import pandas as pd
from snakemake.checkpoints import Checkpoint, Checkpoints
from snakemake.io import expand


def construct_export_files(
    wildcards, manifest: pd.DataFrame, checkpoints: Checkpoints, suffix: str
) -> list:
    """
    Use checkpoint output of linker generation rule to
    determine what files need to be constructed for a
    data export
    """
    res = []
    linker_fn = checkpoints.generate_linker.get().output[0]
    linker_df = pd.read_csv(linker_fn, sep="\t")
    subjectids = manifest.loc[manifest["projectid"] == wildcards.projectid, "sampleid"].to_list()
    targets = linker_df.loc[
        (linker_df["ru"] == wildcards.projectid) & [x in subjectids for x in linker_df["sq"]],
        "output",
    ]
    res = expand(
        "results/export/{projectid}/{file_prefix}.{file_suffix}",
        projectid=wildcards.projectid,
        file_prefix=targets.to_list(),
        file_suffix=suffix,
    )
    return res


def construct_nonexport_files(
    wildcards, manifest: pd.DataFrame, checkpoints: Checkpoints, suffix: str
) -> list:
    """
    Use checkpoint output of linker generation rule to
    determine what files need to be constructed for a
    data (non)export; that is, who isn't mapped to a
    valid external ID
    """
    res = []
    linker_fn = checkpoints.generate_linker.get().output[0]
    linker_df = pd.read_csv(linker_fn, sep="\t")
    subjectids = manifest.loc[manifest["projectid"] == wildcards.projectid, "sampleid"].to_list()
    targets = linker_df.loc[
        (linker_df["ru"] == wildcards.projectid) & [x in subjectids for x in linker_df["sq"]],
        "output",
    ]
    nonexport_targets = []
    for subjectid in subjectids:
        if subjectid not in targets.to_list():
            nonexport_targets.append(subjectid)
    res = expand(
        "results/nonexport/{projectid}/{file_prefix}.{file_suffix}",
        projectid=wildcards.projectid,
        file_prefix=nonexport_targets,
        file_suffix=suffix,
    )
    return res
