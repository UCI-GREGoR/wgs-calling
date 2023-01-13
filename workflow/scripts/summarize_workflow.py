import configparser
import os
import subprocess

import yaml


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


def describe_read_qc(config: dict) -> list:
    """
    Describe the process of ingesting and preprocessing reads,
    up to but not including alignment
    """
    all_lines = ["## Input Read Quality Control"]
    fastqc_version = get_conda_version_string("fastqc")
    fastp_version = get_conda_version_string("fastp")
    ## Describe fastq ingestion. This is a bit verbose, but we're trying to be thorough
    if config["behaviors"]["symlink-fastqs"]:
        res = "Input fastqs were symlinked from external sources hosted locally."
    else:
        res = (
            "Input fastqs were copied from external sources to results/fastqs. Though "
            "this is costly in terms of space utilization, it prevents instances of upstream "
            "files being removed before analysis is complete."
        )
    ## Describe fastqc
    res = (
        res
        + " Fastqs were processed by [fastqc](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) "
        + fastqc_version
        + " with default parameters."
    )
    ## Describe fastp
    res = (
        res
        + " Fastqs were concurrently processed by [fastp](https://github.com/OpenGene/fastp) "
        + fastp_version
        + " with the following "
        "relevant parameter settings: `-q 10 --trim_poly_g --overrepresentation_analysis "
        "--overrepresentation_sampling 100 -5 -3`. Paired, passing reads are passed through fastqc as above "
        "for inspection (as post-trimmed reads) in [multiqc](https://multiqc.info/)."
    )
    all_lines.append(res)
    return all_lines


def describe_alignment(config: dict) -> list:
    """
    Describe the alignment process and its QC, up to but not
    including variant calling
    """
    gatk_version = get_conda_version_string("gatk4")
    alignstats_version = get_conda_version_string("alignstats")
    verify_version = get_conda_version_string("verifybamid2")
    somalier_version = get_conda_version_string("somalier")
    all_lines = ["## Read Alignment"]
    res = "Reads are aligned to human reference genome {} with {}.".format(
        config["genome-build"], config["behaviors"]["aligner"]
    )
    res = (
        res
        + " Aside from standard defaults, softclipping of alts (-Y) and a fixed -K are provided to stabilize output across "
        "repeated runs."
    )
    if config["behaviors"]["trim-adapters-before-alignment"] is True:
        res = (
            res
            + " Following best practices, reads trimmed for the residual Illumina adapter content are used for "
            "alignment. Note that this differs from prior CLO workflow behavior but matches field standards."
        )
    else:
        res = (
            res
            + " Untrimmed reads (that is, not the output of fastp) are "
            "used for alignment. This does not match field standards, but effective soft clipping will likely limit "
            "the impact of this discrepancy on downstream performance."
        )
    res = (
        res
        + " Following best practices, duplicate reads are annotated with [GATK](https://github.com/broadinstitute/gatk) "
        + gatk_version
        + " MarkDuplicates. "
        " Furthermore, duplicate reads "
        "are actually removed from the bam output of this step. This improves compatibility with downstream tools that cannot handle "
        "marked duplicate annotations (e.g. [Sentieon products](https://support.sentieon.com/versions/201808.06/manual/DNAscope_usage/dnascope/))."
    )
    res = (
        res
        + " Again following best practices, quality scores are recalibrated with GATK "
        + gatk_version
        + " BaseRecalibrator/ApplyBQSR."
        "This step was not present in the prior CLO workflow."
    )
    all_lines.append(res)
    res = (
        "Following BQSR, aligned reads are analyzed with GATK "
        + gatk_version
        + " CollectMultipleMetrics, modules "
        "CollectAlignmentSummaryMetrics, CollectInsertSizeMetrics, QualityScoreDistribution, CollectSequencingArtifactMetrics, "
        "CollectQualityYieldMetrics; CollectGcBiasMetrics; and CollectWgsMetrics."
    )
    if config["behaviors"]["trim-adapters-before-alignment"] is True:
        res = (
            res
            + " CollectMultipleMetrics is run with INCLUDE_UNPAIRED active, but "
            "due to the use of fastp before alignment, these reads should be prefiltered."
        )
    else:
        res = (
            res
            + " CollectMultipleMetrics is run with INCLUDE_UNPAIRED active."
        )
    res = (
        res
        + " CollectWgsMetrics is run genome-wide. In the legacy implementation of the pipeline, this step was run "
        "exclusively in the targeted intervals of the NA12878 'reportable regions.' Note that this selection "
        "of regions based on a single HG reference subject is neither appropriate nor evidently originally intended for prior CLO workflow. "
        "This functionality is flagged for additional consideration at a later date."
    )
    res = (
        res
        + " Reads are further analyzed with [alignstats](https://github.com/jfarek/alignstats) "
        + alignstats_version
        + "; [VerifyBamID](https://github.com/Griffan/VerifyBamID) "
        + verify_version
        + "; and [somalier](https://github.com/brentp/somalier) "
        + somalier_version
        + ". Note that currently self-reported sex information is pulled from the upstream 'labbook' sheet, and as such "
        "the completeness of this sexcheck assessment is only as good as the completeness of the 'labbook' annotations. "
        "Furthermore, the relatedness information encoded in the subject ID is not parsed out for direct relatedness checking, "
        "as this is flagged for inclusion in a downstream inter-flowcell QC workflow."
    )
    all_lines.append(res)
    return all_lines


def describe_snv_calling(config: dict) -> list:
    """
    Describe the calling process for SNVs from bams
    """
    bcftools_version = get_conda_version_string("bcftools", True)
    all_lines = ["## SNV Calling"]
    res = ""
    if config["behaviors"]["snv-caller"] == "octopus":
        res = (
            "SNVs were called from aligned bams using octopus. Due to the orphaning of the octopus project, the version of octopus run "
            "is built from [this repository's](https://github.com/invitae-internal/invitae-octopus) commit 17a597d192bcd5192689bf38c5836a98b824867a. "
            "This does not represent an official release version of the software, but rather represents the most recent representation of "
            "the code base before the abandonment of the project. Confidence in this version of the software should be limited. Note that "
            "prior runs of the CLO workflow were using an unreleased version of octopus, but a different unreleased version that the one used "
            "in this workflow; as such, an undefined set of changes may be introduced into the output based on this discrepancy."
        )
    elif config["behaviors"]["snv-caller"] == "deepvariant":
        res = (
            "SNVs were called from aligned bams using [deepvariant](https://github.com/google/deepvariant), with the authors' docker image version "
            + str(config["parameters"]["deepvariant"]["docker-version"])
            + ". Runs were split into "
            + str(config["parameters"]["deepvariant"]["number-shards"])
            + " shards based on the authors' recommendations. To improve integration with UGE/HPC embarrassingly parallel models, the deepvariant "
            "workflow was run as successive steps `make_examples`, `call_variants`, and `postprocess_variants`, followed by a concatenation across shards "
            "using [bcftools](https://github.com/samtools/bcftools) " + bcftools_version + "."
        )
    all_lines.append(res)
    return all_lines


def describe_manta(config: dict, manta_config: str) -> str:
    """
    Describe run configuration for manta
    """
    manta_version = get_conda_version_string("manta")
    all_lines = ["### Manta"]
    cp = configparser.ConfigParser()
    cp.read(manta_config)
    if "manta" in cp:
        res = (
            "[Manta](https://github.com/Illumina/manta) "
            + manta_version
            + " was run with default settings, with the following exceptions."
        )
        all_lines.append(res)
        for setting in cp["manta"].keys():
            if setting == "enableremotereadretrievalforinsertionsingermlinecallingmodes":
                if cp["manta"][setting] == "0":
                    res = (
                        "  - Manual observation has shown that, in spite of internal embarrassingly parallel dispatch, "
                        "Manta gets stuck waiting for hours for a single process to return. "
                        "A [GitHub issue](https://github.com/Illumina/manta/issues/130) concerning this problem suggested "
                        "disabling the parameter 'enableRemoteReadRetrievalForInsertionsInGermlineCallingModes', which "
                        "evidently may degrade sensitivity of certain types of insertions, but by observation drastically "
                        "improves runtime (e.g. reduces runtime by ~75%)."
                    )
                else:
                    res = (
                        "  - The current manta configuration has 'enableRemoteReadRetrievalForInsertionsInGermlineCallingModes' "
                        "activated. This is fine, though observation and a "
                        "[GitHub issue](https://github.com/Illumina/manta/issues/130) suggest that disabling this setting can "
                        "drastically improve overall runtime at the cost of some degradation of sensitivity of certain types "
                        "of insertions."
                    )
            else:
                res = "  - " + setting + ": " + cp["manta"][setting]
            all_lines.append(res)
    else:
        res = (
            "[Manta](https://github.com/Illumina/manta) "
            + manta_version
            + " was run with default settings."
        )
        all_lines.append(res)
    return all_lines


def describe_tiddit(config: dict) -> str:
    """
    Describe run configuration for tiddit
    """
    tiddit_version = get_conda_version_string("tiddit")
    all_lines = ["### TIDDIT"]
    if tiddit_version.startswith("3"):
        res = (
            "[TIDDIT](https://github.com/SciLifeLab/TIDDIT) "
            + tiddit_version
            + " was run with default settings. "
            " "
            "There are certain incompatibilities between TIDDIT >=3 TIDDIT 2. The relevant "
            "workflow modifications are as follows:"
        )
        all_lines.append(res)
        res = (
            "  - TIDDIT uses an internal python fasta parser that implicitly requires the only instance of the '>' character to be "
            "the first character of the sequence description line. "
            "Conservatively, all fastas are sanitized to remove extra '>' characters when ingested into the pipeline."
        )
        all_lines.append(res)
        res = (
            "  - TIDDIT features local realignment of reads, and for this purpose calls bwa. "
            "By default,, the reference fasta has indices that are expected by bwa-mem2, specifically the "
            "modified `*.bwt.2bit.64` file, which differs from the older `*.bwt` file expected by traditional bwa. When the `*.bwt` file is absent, "
            "TIDDIT >=3 fails with a cryptic error. Based on this observation, the `*.bwt` file is separately prepared whenever TIDDIT >=3 is "
            "invoked in this version of the workflow."
        )
        all_lines.append(res)
    else:
        res = "TIDDIT " + tiddit_version + " was run with default settings."
        all_lines.append(res)
    return all_lines


def describe_svaba(config: dict) -> str:
    """
    Describe run configuration for svaba
    """
    svaba_version = get_conda_version_string("svaba")
    all_lines = ["### SvABA"]
    res = (
        "[SvABA](https://github.com/walaj/svaba) "
        + svaba_version
        + " was run with default settings."
    )
    all_lines.append(res)
    return all_lines


def describe_lumpy(config: dict) -> str:
    """
    Describe run configuration for lumpy
    """
    lumpy_version = get_conda_version_string("lumpy-sv")
    smoove_version = get_conda_version_string("smoove")
    all_lines = ["### lumpy"]
    res = (
        "[lumpy](https://github.com/arq5x/lumpy-sv) "
        + lumpy_version
        + " was run according to its project's instructions, under the wrapper program "
        "[smoove](https://github.com/brentp/smoove) "
        + smoove_version
        + ". As this pipeline is focused on single sample calling, only the first "
        "`smoove call` step is run for each sample. The relevant exclusion bedfile for this reference genome is "
        "[here](" + config["lumpy"][config["genome-build"]]["exclude-bed"] + ")."
    )
    all_lines.append(res)
    return all_lines


def describe_delly(config: dict) -> str:
    """
    Describe run configuration for delly
    """
    delly_version = get_conda_version_string("delly")
    all_lines = ["### Delly"]
    res = (
        "[Delly](https://github.com/dellytools/delly) "
        + delly_version
        + " was run with default settings."
    )
    all_lines.append(res)
    return all_lines


def describe_sv_calling(config: dict, manta_config: str) -> list:
    """
    Describe the calling process for SVs from bams
    """
    descriptions = {}
    descriptions["manta"] = describe_manta(config, manta_config)
    descriptions["tiddit"] = describe_tiddit(config)
    descriptions["svaba"] = describe_svaba(config)
    descriptions["lumpy"] = describe_lumpy(config)
    descriptions["delly"] = describe_delly(config)
    duphold_version = get_conda_version_string("duphold")
    svdb_version = get_conda_version_string("svdb")
    ensemble_filtering_criteria = []
    ensemble_min_count = config["behaviors"]["sv-ensemble"]["min-count"]
    if ensemble_min_count > 1:
        ensemble_filtering_criteria.append(
            "present in at least {} tools' post-merge output".format(ensemble_min_count)
        )
    if "required-callers" in config["behaviors"]["sv-ensemble"]:
        ensemble_filtering_criteria.append(
            "must be present in {}".format(
                ",".join(config["behaviors"]["sv-ensemble"]["required-callers"])
            )
        )
    if len(ensemble_filtering_criteria) == 0:
        ensemble_filtering_criteria.append("no additional filters requested in user configuration")
    all_lines = ["## SV Calling"]
    res = (
        "This pipeline is designed to perform ensemble SV calling, combining and filtering the "
        "results of multiple short read callers. This instance of the pipeline is configured to "
        "run: {" + ", ".join(config["behaviors"]["sv-callers"]) + "}."
    )
    if len(config["behaviors"]["sv-callers"]) == 1:
        res = (
            res
            + " Note that because only one caller was used in this run, the benefits of ensemble "
            "calling will be very limited, and the results will reflect the single tool's results with "
            "duphold filters applied as appropriate."
        )
    all_lines.append(res)
    for caller in config["behaviors"]["sv-callers"]:
        all_lines.extend(descriptions[caller])
    all_lines.append("### Depth Annotation and Filtering")
    res = (
        "The results of all SV callers were annotated with depth annotations from [duphold](https://github.com/brentp/duphold) "
        + duphold_version
        + ". Per the tool's recommendations, variants were filtered by the following "
        "criteria: deletion variants required to have at least 0.7X fold change of variant depth "
        "relative to flanking regions, and duplication variants required to have at least 1.3X fold change "
        "of variant depth relative to bins in the genome with similar GC content. Variants are furthermore "
        "required to be annotated with PASS filter status for tools that emit such data. Variants not fitting "
        "in the above descriptions (e.g. interchromosomal breakends) are preserved at this step."
    )
    all_lines.append(res)
    if len(config["behaviors"]["sv-callers"]) == 1:
        res = (
            "The duphold-filtered results of the tool are then merged within-tool with [svdb](https://github.com/J35P312/SVDB) "
            + svdb_version
            + " with default parameters."
        )
    else:
        res = (
            "The duphold-filtered results of each SV tool are then merged with [svdb](https://github.com/J35P312/SVDB) "
            + svdb_version
            + " with default parameters. In brief, svdb first merges similar variants within each tool's individual "
            "output, and then performs a similar merge of variants across the remaining variants from each tool. After "
            "the merge is complete, this instance of the pipeline is configured to filter variants based on the following "
            "criteria: {" + "; ".join(ensemble_filtering_criteria) + "}."
        )
    all_lines.append("### Variant Ensemble Calling")
    all_lines.append(res)
    return all_lines


def describe_data_release(config: dict) -> None:
    """
    Describe the steps applied to data release
    """
    all_lines = ["## Data Preparation for Release"]
    all_lines.append("### Aligned Reads")
    all_lines.append(
        "Aligned bam files are mapped to generic sample IDs {PMGRCID_LSID_SQID}. For recordkeeping, this pipeline's version is added "
        "to the bam header as a comment (@CO) tag."
    )
    all_lines.append("### SNV Calls")
    all_lines.append(
        "SNV VCFs are mapped to generic sample IDs {PMGRCID_LSID_SQID}, both in filename and in vcf sample header."
        "SNVs are filtered based on the following criteria (derived from [Pedersen _et al._](https://doi.org/10.1038/s41525-021-00227-3):\n"
        "  - filter status PASS\n"
        "  - genotype quality >= 20\n"
        "  - genotype total depth >= 10\n"
        "  - allele balance (heterozygotes) between 0.2 and 0.8, inclusive\n"
        "  - allele balance (homozygous alts) less than 0.04\n"
        "    - this is the only meaningful deviation from the above citation,\n"
        "      and is based on the observation that a meaningful proportion of\n"
        "      calls have exactly one reference read, and at least by observation\n"
        "      it seems like these variants don't deserve to be filtered out at this step\n"
        '  - not intersecting with ENCODE "blacklist V2" regions [here](https://github.com/Boyle-Lab/Blacklist)\n'
        "Multiallelics are split with [bcftools norm -m -both](https://samtools.github.io/bcftools/bcftools.html#norm).\n"
        "For recordkeeping purposes, this pipeline's version is added to the vcf header. For compatibility with Moon,\n"
        "the genome reference code provided in configuration ("
        + config["genome-build"]
        + ") is added to the vcf header."
    )
    return all_lines


def run_summarize_workflow(config: dict, manta_config: str, oname: str) -> None:
    """
    Describe the current workflow structure in prose
    with applicable software version data
    """
    process_description = ["# wgs-pipeline Methods Summary"]
    process_description.append("## Overview")
    process_description.append(
        "This instance of [wgs-pipeline](https://gitlab.com/lightning.auriga1/wgs-pipeline) "
        "is based on git version " + get_pipeline_git_version() + ". What follows is a prose-style "
        "description both of the general methodology of the pipeline, and of the effect of the user "
        "configuration settings currently in effect in this pipeline instance. Software version data "
        "are derived from [conda](https://docs.conda.io/en/latest/) package versions. Note that this "
        "information is only as accurate as the correspondence between the data in this workflow instance "
        "and the settings currently present in the user configuration file (config/config.yaml), and "
        "delinking will occur if these files are retroactively modified after the run is complete."
    )
    ## All pipeline instances run read QC
    process_description.extend(describe_read_qc(config))
    ## Only include alignment description if the user requested it, either directly or implicitly
    if (
        config["behaviors"]["outcome"] == "alignment"
        or config["behaviors"]["outcome"] == "calling"
        or config["behaviors"]["outcome"] == "release"
    ):
        process_description.extend(describe_alignment(config))
    ## Only include calling description if the user requested it
    if config["behaviors"]["outcome"] == "calling" or config["behaviors"]["outcome"] == "release":
        process_description.extend(describe_snv_calling(config))
        process_description.extend(describe_sv_calling(config, manta_config))
    with open(oname, "w") as f:
        f.writelines(["{}\n\n".format(x) for x in process_description])
    ## Only include data release description if the user requested it
    if config["behaviors"]["outcome"] == "release":
        process_description.extend(describe_data_release(config))


run_summarize_workflow(
    snakemake.params["config"],  # noqa: F821
    snakemake.input["manta_config"],  # noqa: F821
    snakemake.output["markdown"],  # noqa: F821
)
