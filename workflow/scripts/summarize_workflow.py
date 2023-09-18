import configparser
import os
import subprocess

import yaml
from jinja2 import Environment, FileSystemLoader


def get_pipeline_git_version() -> str:
    """
    Query the current instance of the workflow to determine the best
    representation of its git version
    """
    commit_id = subprocess.check_output(["git", "rev-parse", "HEAD"], encoding="UTF-8").rstrip()
    commit_branch = subprocess.check_output(
        ["git", "rev-parse", "--abbrev-ref", "HEAD"], encoding="UTF-8"
    ).rstrip()
    try:
        commit_description = subprocess.check_output(
            ["git", "describe", "--tags", commit_branch],
            encoding="UTF-8",
            stderr=subprocess.DEVNULL,
        )
    except subprocess.CalledProcessError:
        commit_description = commit_id
    else:
        commit_description = commit_description.rstrip()
    return commit_description


def get_conda_version_string(pkgname: str, require_first: bool = False):
    """
    Given a conda package name, find out what installed version
    of the package is available in a workflow's conda environments
    """
    conda_dir = ".snakemake/conda"
    yaml_files = [
        "{}/{}".format(conda_dir, y)
        for y in filter(lambda x: x.endswith(".yaml"), os.listdir(conda_dir))
    ]
    candidate_yamls = []
    candidate_tags = []
    for yaml_file in yaml_files:
        with open(yaml_file, "r") as f:
            config = yaml.safe_load(f)
        if pkgname in config["dependencies"]:
            if require_first and config["dependencies"][0] != pkgname:
                continue
            env_prefix = yaml_file.removesuffix(".yaml")
            candidate_yamls.append(env_prefix)
    for env_prefix in candidate_yamls:
        targets = [
            "{}/conda-meta/{}".format(env_prefix, y)
            for y in filter(
                lambda x: x.startswith(pkgname),
                os.listdir("{}/conda-meta".format(env_prefix)),
            )
        ]
        ## let's make the dangerous assumption that the first such result is the one we want.
        ## there are instances of packages with identical prefix strings (e.g. lumpy-sv- and
        ## lumpy-sv-minimal-). in that particular instance, the actual version number is the same;
        ## and while in such instances we'd often expect the version numbers to be the same,
        ## there's no easy guarantee of that
        pkg_tag = (
            os.path.basename(targets[0])
            .removeprefix(pkgname + "-")
            .removesuffix(".json")
            .split("-")[0]
        )
        candidate_tags.append(pkg_tag)
    ## this script does not currently have a method of distinguishing between multiple installations of the same
    ## conda environment that have been installed in the same place due to e.g. development changes to conda env specs.
    ## for the moment, we can at least detect if this has occurred, and not return the entirely wrong answer.
    if len(set(candidate_tags)) == 1:
        return candidate_tags[0]
    elif len(candidate_tags) > 0:
        return "{multiple conflicting versions found in installed conda environments}"
    else:
        return "{not found in installed conda environments}"


def get_all_sv_callers(config: dict) -> list:
    """
    Across all ensemble calling paradigms, get all requested SV callers
    """
    res = []
    for endpoint in config["behaviors"]["sv-endpoints"].keys():
        res.extend(config["behaviors"]["sv-endpoints"][endpoint]["sv-callers"])
    return res


def describe_ensemble_calling_settings(config: dict) -> str:
    """
    Parse config settings into human-readable description of
    SV ensemble calling settings
    """
    ensemble_filtering_criteria = ""
    for endpoint in config["behaviors"]["sv-endpoints"].keys():
        ensemble_filtering_criterion = '\n\nendpoint "' + endpoint + '": '
        sv_callers = config["behaviors"]["sv-endpoints"][endpoint]["sv-callers"]
        ensemble_filtering_criterion += "SV callers: " + ", ".join(sv_callers) + "; "
        ensemble_min_count = config["behaviors"]["sv-endpoints"][endpoint]["sv-ensemble"][
            "min-count"
        ]
        if ensemble_min_count > 1:
            ensemble_filtering_criterion += (
                "present in at least {} tools' post-merge output".format(ensemble_min_count) + "; "
            )
        if "required-callers" in config["behaviors"]["sv-endpoints"][endpoint]["sv-ensemble"]:
            ensemble_filtering_criterion += (
                "must be present in {}".format(
                    ",".join(config["behaviors"]["sv-ensemble"]["required-callers"])
                )
                + "; "
            )
        if config["behaviors"]["sv-endpoints"][endpoint]["sv-remove-breakends"]:
            ensemble_filtering_criterion += "breakends removed from processed results; "
        ensemble_filtering_criteria += ensemble_filtering_criterion
    return ensemble_filtering_criteria.rstrip("; ")


def run_summarize_workflow(template_name, config, manta_config_fn, outfn):
    """
    Initialize variables and pass them into jinja along with template
    """
    env = Environment(loader=FileSystemLoader("."), extensions=["jinja_markdown.MarkdownExtension"])
    jinja_template = env.get_template(template_name)
    function_dict = {"get_conda_version_string": get_conda_version_string}
    cp = configparser.ConfigParser()
    cp.read(manta_config_fn)
    ensemble_filtering_criteria = describe_ensemble_calling_settings(config)
    jinja_template.globals.update(function_dict)
    template_string = jinja_template.render(
        config=config,
        manta_settings=cp,
        ensemble_filtering_criteria=ensemble_filtering_criteria,
        pipeline_git_version=get_pipeline_git_version(),
    )
    with open(outfn, "w") as f:
        f.write(template_string)


run_summarize_workflow(
    snakemake.input["jinja_template"],  # noqa: F821,
    snakemake.params["config"],  # noqa: F821
    snakemake.input["manta_config"],  # noqa: F821
    snakemake.output["html"],  # noqa: F821
)
