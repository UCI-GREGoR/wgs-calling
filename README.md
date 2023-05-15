# Snakemake workflow: WGS Pipeline

This workflow is intended to be the R&D space for the PMGRC WGS analysis pipeline. It will first support existing functionality from variously [marigold](https://github.com/invitae-internal/marigold-pipes), [nextflow-pipelines](https://github.com/invitae-internal/nextflow-pipelines), and the descendants of `nextflow-pipelines` after splitting out into individual repositories. When that back compatibility is complete, additional features will be sandboxed and tested here.

New global targets should be added in `workflow/Snakefile`. Content in `workflow/Snakefile` and the snakefiles in `workflow/rules` should be specifically _rules_; python infrastructure should be composed as subroutines under `lib/` and constructed in such a manner as to be testable with [pytest](https://docs.pytest.org/en/7.2.x/). Rules can call embedded scripts (in python or R/Rmd) from `workflow/scripts`; again, these should be constructed to be testable with pytest or [testthat](https://testthat.r-lib.org/).

## Authors

* Lightning Auriga (@lightning.auriga)

## Usage

### Step 1: Obtain a copy of this workflow

1. Clone this repository to your local system, into the place where you want to perform the data analysis.
```
    git clone git@gitlab.com:lightning.auriga1/wgs-pipeline.git
```

Note that this requires local git ssh key configuration; see [here](https://docs.gitlab.com/ee/user/ssh.html) for instructions as required.

### Step 2: Configure workflow

Configure the workflow according to your needs via editing the files in the `config/` folder. Adjust `config.yaml` to configure the workflow execution, and `manifest.tsv` to specify your sample setup.

The following settings are recognized in `config/config.yaml`. Note that each reference data option below exists under an arbitrary tag denoting desired reference genome build. This tag is completely arbitrary and will be used to recognize the requested build for the current pipeline run.

|Configuration Setting|Description|
|---|---|
|`manifest`|relative path to run manifest|
|`sample-logbook`|local Excel spreadsheet clone of sample manifest information from Google docs.<br>This is upstream input. a local cloned file is preferred due to the possibility of uncontrolled upstream changes|
|`sample-linking`|two column tsv linker files with sample annotation data|
||`sex`: path to linker containing internal subject ID and self-reported sex data|
||`external-ids`: path to linker containing internal subject ID and external IDs|
|`multiqc-read-config`|relative path to configuration settings for pre-alignment multiQC report|
|`multiqc-alignment-config`|relative path to configuration settings for post-alignment multiQC report|
|`genome-build`|requested genome reference build to use for this analysis run. this should match the tags used in the reference data blocks below.|

The following settings are nested under the key `behaviors` and are user-configurable modifiers to how the pipeline will run.

|Behavior|Description|
|---|---|
|`use-containers`|whether to, when possible, use either the docker or singularity image for each rule, or instead the rule-specific conda environment. See discussion below for how to choose this setting, and how it interacts with snakemake invocations.
|`aligner`|which alignment tool to use. permitted values: `bwa`, `bwa-mem2`|
|`snv-caller`|which calling tool to use for SNVs. permitted values: `deepvariant`|
|`sv-callers`|which calling tool(s) to use for SVs. at least one should be specified. permitted values: `manta`, `tiddit`, `svaba`, `delly`, `lumpy`
|`sv-ensemble`|settings controlling SV ensemble calling. note that the below settings can be applied in combination|
||`min-count`: the minimum number of tools' outputs in which a variant (or something similar nearby) must appear to survive ensemble filtering|
||`required-callers`: a list of which tools, if any, a variant absolutely must appear in to survive ensemble filtering
|`sv-remove-breakends`|whether or not to filter `SVTYPE=BND` variants from ensemble calling output
|`outcome`|which endpoint to run to. permitted values: `fastqc` (for read QC only); `alignment`; or `calling`; or `release` to prepare results for distribution|
|`symlink-fastqs`|whether to copy (no) or symlink (yes) input fastqs into workspace. symlinking is faster and more memory-efficient, but less reproducible, as the upstream files may vanish leaving no way to regenerate your analysis from scratch.|
|`trim-adapters-before-alignment`|whether to use adapter trimmed fastq output of `fastp` as input to aligner. permitted values: `yes`, `no`, or `legacy`. legacy behavior for this option is to not use trimmed output for alignment.|
|`assume-last-sample-sex`|upstream convention is to include a low-depth NA24385 as the final sample in every flowcell. that sample is not annotated in the logbook in a way that the parser understands. to maintain compatibility with Somalier sexcheck, set this to `female`. if no sex should be assumed for the sample, delete this configuration option|
|`export-directory`|top-level path to where output files should be moved after release run is complete. delete this option to disable.
|`export-s3`|parameters for controlling optional upload to s3|
||`bucket-name`: name of s3 bucket to which to sync data|
||`profile-name`: optional name of aws profile to use for data sync|

The following tool-specific parameters are nested under the key `parameters`.

|Tool|Recognized Parameters|
|---|---|
|`bwa`|parameters specific to [bwa](https://bio-bwa.sourceforge.net/)|
||`K`: chunk size parameter. `bwa` defaults this to `1e7*{threads}`, but to maintain consistency independent of thread count, this is manually fixed. higher numbers improve runtime at the cost of (marginally) increased RAM usage|
|`bwa-mem2`|parameters specific to [bwa-mem2](https://github.com/bwa-mem2/bwa-mem2)|
||`K`: chunk size parameter. see the corresponding `bwa` option for details|
|`deepvariant`|parameters specific to [deepvariant](https://github.com/google/deepvariant)|
||`number-shards`: how many shards to break calling into. needs to be at most the number of available threads in the submission queue|
||`docker-version`: which docker tag to use when pulling the official DeepVariant docker image|
|`manta`|parameters specific to [Manta](https://github.com/Illumina/manta)|
||`config-ini`: relative path to Manta local configuration data. see Manta documentation for more specific information about permitted values in this config file. the exposure of this file, which is passed to `configManta.py`, is done in anticipation of toggling settings that can reduce Manta's substantial runtime, e.g. `enableRemoteReadRetrievalForInsertionsInGermlineCallingModes`
|`tiddit`|parameters specific to [TIDDIT](https://github.com/SciLifeLab/TIDDIT)|
||`min-contig-size`: minimum size of contigs on which to make calls, in bases. a minimum size of 2000000 will remove the remaining non-standard contigs present in the no-alt version of GRCh38, for example. unfortunately, TIDDIT does not directly support calling regions|

The following general genome reference files are nested under the key `parameters`, and are expected to be available to multiple tools in the pipeline.

|Annotation Type|Description|
|---|---|
|`fasta`|human sequence fasta file.|
|`exclusion-regions-bed`|set of bed regions to be excluded from output SNV data. expected to be from https://doi.org/10.1038/s41598-019-45839-z|
|`exons-gtf-gz`|gene annotations for genome build. used to extract a bedfile of exon positions|

The following reference files are split out by individual tool.

|Tool|Annotation Type|Description|
|---|---|---|
|`bsqr`||reference files for base quality score recalibration from GATK4|
||`known-indels-vcf-gz`|VCF of "gold standard" indels for BQSR. intended to be pulled from Broad's cloud files|
||`known-indels-vcf-gz-tbi`|tabix index for above known indels vcf|
||`dbsnp138-vcf`|VCF of dbSNP variation for BQSR. intended to be pulled from Broad's cloud files|
||`dbsnp138-vcf-idx`|.idx index file for above dbSNP vcf|
|`collectwgsmetrics`||reference files specific for GATK4 post-alignment QC utility `CollectWgsMetrics`|
||`reportable-regions`|set of bed intervals to be considered for `CollectWgsMetrics` output. note that the legacy Nextflow workflow unconditionally applied the NA12878 reportable regions set to all samples. this is likely in error, and as such will need to be changed|
|`verifybamid2`||reference data files specific to [VerifyBamID2](https://github.com/Griffan/VerifyBamID)|
||`db-V`|filename for assorted Verify annotation files|
||`db-UD`|filename for assorted Verify annotation files|
||`db-mu`|filename for assorted Verify annotation files|
||`db-bed`|filename for assorted Verify annotation files|
|`deepvariant`||reference data files specific to [deepvariant](https://github.com/google/deepvariant)|
||`calling-ranges`|list of files containing chromosomal intervals for embarrassingly parallel analysis. these are currently effectively placeholders that just contain the standard human chromosomes, but at some point, this will be extended to contain actual balanced calling intervals|
|`lumpy`||reference data files specific to [lumpy](https://github.com/arq5x/lumpy-sv) or [smoove](https://github.com/brentp/smoove)|
||`exclude-bed`|set of hard-to-call bed intervals to exclude from analysis. these files are intended to be the ones listed in the tools' documentation, but can be customized if desired|
|`delly`||reference data files specific to [delly](https://github.com/dellytools/delly)|
||`exclude-bed`|set of hard-to-call bed intervals to exclude from analysis. these files are intended to be the ones listed in the tools' documentation, but can be customized if desired|
|`svaba`||reference data files specific to [svaba](https://github.com/walaj/svaba)|
||`exclude-bed`|set of hard-to-call bed intervals to exclude from analysis. though this option is exposed by SvABA, there is no recommended file to use here, so one might just use the lumpy files, for example|
|`manta`||reference data files specific to [manta](https://github.com/Illumina/manta)|
||`calling-range-bed-gz`|bgzip-compressed bedfile of valid calling ranges, e.g. autosomes. manta recommends not providing very many regions here, as it evidently causes problems with its task dispatch heuristics
||`calling-range-bed-gz-tbi`|tabix index of above compressed bedfile|

The following columns are expected in the run manifest, by default at `config/manifest.tsv`:

|Manifest column|Description|
|---|---|
|`projectid`|run ID, or other desired grouping of sequencing samples. this will be a subdirectory under individual tools in `results/`|
|`sampleid`|sequencing ID for sample|
|`r1`|R1 fastq.gz file for sample|
|`r2`|R2 fastq.gz file for sample|
|`lane`|(optional) sequencing lane code, with `L00` prefix. if not specified, will be assumed to be `L001`|

### Step 3: Install Snakemake

Install Snakemake using [conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html):

    conda create -c bioconda -c conda-forge -n snakemake snakemake

For installation details, see the [instructions in the Snakemake documentation](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html).

### Step 4: Execute workflow

Activate the conda environment:

    conda activate snakemake

Test your configuration by performing a dry-run via

    snakemake --use-conda -n

Execute the workflow locally via

    snakemake --use-conda --cores $N

using `$N` cores or run it in a cluster environment via

    snakemake --use-conda --use-singularity --profile sge-profile --cluster-config config/cluster.yaml --jobs 100

See the [Snakemake documentation](https://snakemake.readthedocs.io/en/stable/executable.html) for further details.

#### Cluster profiles

Snakemake interfaces with job schedulers via _cluster profiles_. For running jobs on SGE, you can use
the cookiecutter template [here](https://github.com/Snakemake-Profiles/sge).



#### How to choose the pipeline's rule-specific dependency behavior

This pipeline is designed to be run using either conda or docker/singularity/apptainer to manage rule-specific dependencies.
Snakemake's usual method for selecting between these options is to either specify `--use-conda` or `--use-singularity`
during the `snakemake` command line invocation. However, due to a lack of a functional conda environment for DeepVariant,
pure conda mode is not possible when using DeepVariant for variant calling.

To get around this discrepancy, this workflow should always be invoked with `snakemake --use-conda --use-singularity`.
The user can then further control whether they want pure-container or (almost) pure-conda mode by setting the
userspace configuration `use-containers` in `config/config.yaml`. In order to use containers, a local copy of the
apptainer images built for this workflow is required; the method of acquiring these containers is TBD but will
probably involve an S3 pull.

### Step 5: Investigate results

#### Read Quality Control

Quality control data from fastqc and fastp are aggregated in a multiqc report at `results/multiqc/{projectid}/multiqc.fastq.html`.
This version of the quality control report is split by lane within the report, if per-lane fastqs have been provided
and annotated in the manifest.

#### Alignment

Quality control data from the above read QC tools as well as somalier, verifybamid, alignstats,  and assorted gatk4/picard analysis tools
are aggregated in a multiqc report at `results/multiqc/{projectid}/multiqc.alignment.html`. This version of the quality report is currently
under active modification, and should only be used for considering alignment QC results (as opposed to pre-alignment read QC) until further notice.

#### Variant Calling

All variant calling in this pipeline is conducted per-subject, ignoring batch data. SNV calls from the user-configured tool
(e.g. DeepVariant) are aggregated in `results/{toolname}/{projectid}/{sampleid}.sorted.vcf.gz`. SV calls from
ensemble calling based on user-configured tools and exclusion criteria are aggregated in `results/final/{projectid}/{sampleid}.sv.vcf.gz`.
Note that these paths and filenames are subject to change before stabilization of the workflow.

#### GVCFs for Batch Calling

If the user has selected a compatible tool for SNV calling, gvcf files per-subject will be collected at
`results/deepvariant/{projectid}/{sampleid}.g.vcf.gz`. These files are not used in the pipeline itself,
and represent the raw output of `deepvariant postprocess_variants`. At some point, these will likely receive
further processing in anticipation of use with e.g. GLnexus in a different pipeline. Also at some point,
I'll probably add a flag for disabling the production of gvcf output, but that's not urgent.

#### Optional: emit methods description and software version data

Some users may be interested in a specific breakdown of workflow methods, relevant software versions pulled from
conda, and the effects of certain important user configuration settings. This information can be generated upon
completion of a workflow with the command `snakemake -j1 --use-conda results/reports/methods_summary.html`.
This will create an additional output file, `results/reports/methods_summary.html`, that contains the best
description of the workflow as was actually run. Note that this information focuses on methodology as opposed to
results, and only requires the relevant conda environments be installed; so if you want to predict what the
workflow will do before actually running it, complete user configuration, run
`snakemake -j1 --use-conda --conda-create-envs-only`, and then run the aforementioned snakemake command
target to generate a markdown description of what _would_ happen if the pipeline were deployed.

#### Data Release

When the pipeline is run in `release` mode, postprocessed output files for each flowcell will be emitted under
`results/export/{flowcell_id}`. The following files will be present, with modifications as annotated:

- aligned reads, represented as lossless crams
  - the source reference fasta used for these files is both in the cram header as a `@CO` and also
    annotated in the output methods html
- `crai` index files for the above cram files
- SNV vcfs, bgzip-compressed, with the following modifications (derived from [Pedersen _et al._](https://doi.org/10.1038/s41525-021-00227-3):
  - only FILTER=PASS variants
  - multiallelics split to biallelics
  - GQ >= 20
  - DP >= 10
  - for heterozygotes, allele balance on `[0.2, 0.8]`
  - for homozygous alts, allele balance `<= 0.04`
  - variants intersected with configured exclusion regions removed
- tabix indices for the above vcfs
- md5 sums for all above files
- a plaintext `manifest.tsv` containing a list of the above data files
- an immutable `methods_summary.html` representing the rendered version of the above methods description for the version of the pipeline that created the release
- SVs TBD

#### (optional) Move Data to External Release Directory

A convenience rule exists for exporting data out of the workflow and into some sort of persistent storage location.
This rule is executed by running `snakemake -j1 export_data`. This rule has the following assumptions:

- it assumes the current working directory has a JIRA ticket formatted `RT-\d{4}` somewhere in its absolute path
- it assumes the current working directory has a flowcell ID formatted `RU\d{5}` somewhere in its absolute path
- it assumes a target export directory has been defined in `config/config.yaml` under `behaviors::export-directory`
- it assumes the workflow has already been run to completion in `release` mode

If the above assumptions are met, the rule will do the following:

- construct a subdirectory under `export-directory` named `{JIRA ticket}/{flowcell ID}`
- move all `*vcf.gz*` and `*cram*` files from `results/export` to that directory
- move the methods summary from `results/export` to that directory
- edit the checksum files that were in `results/export` to no longer contain the relative paths specific to the workflow results directory structure
- run `md5sum -c` on all files with checksums (cram, crai, vcf.gz, tbi) and report the results to `results/export/md5_checks.txt`

These behaviors are controlled by the utility shell script in `workflow/scripts/export_data.bash`, which can be called outside of the pipeline if
desired as `./export_data.bash {export_directory} {md5_check_outfile.txt}`. This is just included for convenience as this file export step is required
for all runs of the pipeline for my particular use case, and can be ignored if desired.

### Step 6: Commit changes

Whenever you change something, don't forget to commit the changes back to your github copy of the repository:

    git commit -a
    git push

### Step 7: Obtain updates from upstream

Whenever you want to synchronize your workflow copy with new developments from upstream, do the following.

1. Once, register the upstream repository in your local copy: `git remote add -f upstream git@gitlab.com:lightning.auriga1/wgs-pipeline.git` or `upstream https://gitlab.com/lightning.auriga1/wgs-pipeline.git` if you do not have setup ssh keys.
2. Update the upstream version: `git fetch upstream`.
3. Create a diff with the current version: `git diff HEAD upstream/default workflow > upstream-changes.diff`.
4. Investigate the changes: `vim upstream-changes.diff`.
5. Apply the modified diff via: `git apply upstream-changes.diff`.
6. Carefully check whether you need to update the config files: `git diff HEAD upstream/default config`. If so, do it manually, and only where necessary, since you would otherwise likely overwrite your settings and samples.


### Step 8: Contribute back

In case you have also changed or added steps, please consider contributing them back to the original repository. This project follows git flow; feature branches off of dev are welcome.

1. [Clone](https://docs.gitlab.com/ee/gitlab-basics/start-using-git.html) the fork to your local system, to a different place than where you ran your analysis.
2. Check out a branch off of dev:
```
git fetch
git checkout dev
git checkout -b your-new-branch
```
3. Make whatever changes best please you to your feature branch.
4. Commit and push your changes to your branch.
5. Create a [merge request](https://docs.gitlab.com/ee/user/project/merge_requests/) against dev.

## Testing

Testing infrastructure for embedded python and R scripts is installed under `lib/` and `workflow/scripts/`. Additional testing
coverage for the Snakemake infrastructure itself should be added once the workflow is more mature ([see here](https://github.com/lightning.auriga/snakemake-unit-tests)).

### Python testing with `pytest`
The testing under `lib/` is currently functional. Partial testing exists for the builtin scripts under `workflow/scripts`: the new utilities
for this implementation are tested, but some code inherited from the legacy pipeline(s) is not yet covered. To run the tests, do the following (from top level):

```bash
mamba install pytest-cov
pytest --cov=lib --cov=workflow/scripts lib workflow/scripts
```


### R testing with `testthat`
The testing under `workflow/scripts` is currently functional. The tests can be run with the utility script `run_tests.R`:

```bash
Rscript ./run_tests.R
```

To execute the above command, the environment must have an instance of R with appropriate libraries:

```bash
mamba install -c bioconda -c conda-forge "r-base>=4" r-testthat r-covr r-r.utils r-desctools
```
