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

- `manifest`: relative path to run manifest
- `multiqc-config`: relative path to configuration settings for post-alignment multiQC report
- `genome-build`: requested genome reference build to use for this analysis run. this should match the tags used in the reference data blocks below.
- `behaviors`: user-configurable modifiers to how the pipeline will run
  - `aligner`: which alignment tool to use. permitted values: `bwa-mem2`
  - `snv-caller`: which calling tool to use for SNVs. permitted values: `octopus`
  - `sv-callers`: which calling tool(s) to use for SVs. at least one should be specified. permitted values: `manta`, `tiddit`, `svaba`, `delly`, `lumpy`
  - `sv-ensemble`: settings controlling SV ensemble calling. note that the below settings can be applied in combination
	- `min-count`: the minimum number of tools' outputs in which a variant (or something similar nearby) must appear to survive ensemble filtering
	- `required-callers`: a list of which tools, if any, a variant absolutely must appear in to survive ensemble filtering
  - `outcome`: which endpoint to run to. permitted values: `fastqc` (for read QC only); `alignment`; or `calling`; or `release` to prepare results for distribution
  - `symlink-fastqs`: whether to copy (no) or symlink (yes) input fastqs into workspace. symlinking is faster and more memory-efficient, but
    less reproducible, as the upstream files may vanish leaving no way to regenerate your analysis from scratch.
  - `trim-adapters-before-alignment`: whether to use adapter trimmed fastq output of `fastp` as input to aligner.
    permitted values: `yes`, `no`, or `legacy`. legacy behavior for this option is to not use trimmed output for alignment.
- `parameters`: tool-specific parameters. note that this section is a work in progress, somewhat more than the rest
  - `deepvariant`: parameters specific to [deepvariant](https://github.com/google/deepvariant)
    - `number-shards`: how many shards to break calling into. needs to be at most the number of available threads in the submission queue
	- `docker-version`: which docker tag to use when pulling the official DeepVariant docker image
  - `manta`: parameters specific to [Manta](https://github.com/Illumina/manta)
    - `config-ini`: relative path to Manta local configuration data. see Manta documentation for more specific information about permitted
	  values in this config file. the exposure of this file, which is passed to `configManta.py`, is done in anticipation of toggling
	  settings that can reduce Manta's substantial runtime, e.g. `enableRemoteReadRetrievalForInsertionsInGermlineCallingModes`
  - `tiddit`: parameters specific to [TIDDIT](https://github.com/SciLifeLab/TIDDIT)
	- `min-contig-size`: minimum size of contigs on which to make calls, in bases. a minimum size of 2000000 will remove the remaining
	  non-standard contigs present in the no-alt version of GRCh38, for example. unfortunately, TIDDIT does not directly support calling regions
- `references`: human genome reference data applicable to multiple tools
  - `fasta`: human sequence fasta file
  - note that the other bwa-style index files attached to this fasta used to be imported by the nextflow workflow. however, presumably by accident,
    these annotation files were getting pulled from various different directories in a way that suggested that they might be delinked from their
	source fasta. in reality, the source reference fastas were probably the same; but to avoid any difficulties downstream, now only the fasta
	itself is pulled in from remote, and the index files are regenerated. this also substantially cleans up the configuration.
  - `exclusion-regions-bed`: set of bed regions to be excluded from output SNV data. expected to be from https://doi.org/10.1038/s41598-019-45839-z
- `bsqr`: reference files for base quality score recalibration from GATK4
  - `known-indels-vcf-gz`: VCF of "gold standard" indels for BQSR. intended to be pulled from Broad's cloud files
  - `known-indels-vcf-gz-tbi`: tabix index for above known indels vcf
  - `dbsnp138-vcf`: VCF of dbSNP variation for BQSR. intended to be pulled from Broad's cloud files
  - `dbsnp138-vcf-idx`: .idx index file for above dbSNP vcf
- `collectwgsmetrics`: reference files specific for GATK4 post-alignment QC utility `CollectWgsMetrics`
  - `reportable-regions`: set of bed intervals to be considered for `CollectWgsMetrics` output. note that the legacy Nextflow workflow
    unconditionally applied the NA12878 reportable regions set to all samples. this is likely in error, and as such will need to be changed
- `dnascope`: reference data files specific to [Sentieon DNAscope](https://support.sentieon.com/manual/DNAscope_usage/dnascope/)
  - `model`: DNAscope model file
  - `dbsnp-vcf-gz`: dbSNP backend vcf.gz file
- `verifybamid2`: reference data files specific to [VerifyBamID2](https://github.com/Griffan/VerifyBamID)
  - `db-V`: filename for assorted Verify annotation files
  - `db-UD`: filename for assorted Verify annotation files
  - `db-mu`: filename for assorted Verify annotation files
  - `db-bed`: filename for assorted Verify annotation files
- `octopus`: reference data files specific to [octopus](https://github.com/luntergroup/octopus)
  - `forest-model`: forest model annotation file for `--forest-model` [note: this is independent of genome build (probably)]
  - `error-model`: error model annotation file for `--sequence-error-model` [note: this is independent of genome build (probably)]
  - `skip-regions`: region annotation for `--skip-regions-file`
  - `calling-ranges`: list of files containing chromosomal intervals for embarrassingly parallel analysis. these are currently effectively
    placeholders that just contain the standard human chromosomes, but at some point, this will be extended to contain actual
	balanced calling intervals
- `deepvariant`: reference data files specific to [deepvariant](https://github.com/google/deepvariant)
  - `calling-ranges`: list of files containing chromosomal intervals for embarrassingly parallel analysis. these are currently effectively
    placeholders that just contain the standard human chromosomes, but at some point, this will be extended to contain actual
	balanced calling intervals
- `lumpy`: reference data files specific to [lumpy](https://github.com/arq5x/lumpy-sv) or [smoove](https://github.com/brentp/smoove)
  - `exclude-bed`: set of hard-to-call bed intervals to exclude from analysis. these files are intended to be the ones listed in the tools' documentation,
    but can be customized if desired
- `delly`: reference data files specific to [delly](https://github.com/dellytools/delly)
  - `exclude-bed`: set of hard-to-call bed intervals to exclude from analysis. these files are intended to be the ones listed in the tools' documentation,
    but can be customized if desired
- `svaba`: reference data files specific to [svaba](https://github.com/walaj/svaba)
  - `exclude-bed`: set of hard-to-call bed intervals to exclude from analysis. though this option is exposed by SvABA, there is no recommended file
    to use here, so one might just use the lumpy files, for example
- `manta`: reference data files specific to [manta](https://github.com/Illumina/manta)
  - `calling-range-bed-gz`: bgzip-compressed bedfile of valid calling ranges, e.g. autosomes. manta recommends not providing very many regions here,
    as it evidently causes problems with its task dispatch heuristics
  - `calling-range-bed-gz-tbi`: tabix index of above compressed bedfile

The following columns are expected in the run manifest, by default at `config/manifest.tsv`:
- `projectid`: run ID, or other desired grouping of sequencing samples. this will be a subdirectory under individual tools in `results/`
- `sampleid`: sequencing ID for sample
- `r1`: R1 fastq.gz file for sample
- `r2`: R2 fastq.gz file for sample
- `lane`: (optional) sequencing lane code, with `L00` prefix. if not specified, will be assumed to be `L001`

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

    snakemake --use-conda --profile sge-profile --cluster-config config/cluster.yaml --jobs 100

See the [Snakemake documentation](https://snakemake.readthedocs.io/en/stable/executable.html) for further details.

#### Cluster profiles

Snakemake interfaces with job schedulers via _cluster profiles_. For running jobs on SGE, you can use
the cookiecutter template [here](https://github.com/Snakemake-Profiles/sge).



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
(e.g. DeepVariant, Octopus) are aggregated in `results/{toolname}/{projectid}/{sampleid}.sorted.vcf.gz`. SV calls from
ensemble calling based on user-configured tools and exclusion criteria are aggregated in `results/final/{projectid}/{sampleid}.sv.vcf.gz`.
Note that these paths and filenames are subject to change before stabilization of the workflow.

#### (DeepVariant only) GVCFs for Batch Calling

If the user has selected DeepVariant for SNV calling, gvcf files per-subject will be collected at
`results/deepvariant/{projectid}/{sampleid}.g.vcf.gz`. These files are not used in the pipeline itself,
and represent the raw output of `deepvariant postprocess_variants`. At some point, these will likely receive
further processing in anticipation of use with e.g. GLnexus in a different pipeline. Also at some point,
I'll probably add a flag for disabling the production of gvcf output, but that's not urgent.

gvcf output is not supported by Octopus and so is not possible in this pipeline.

#### Optional: emit methods description and software version data

Some users may be interested in a specific breakdown of workflow methods, relevant software versions pulled from
conda, and the effects of certain important user configuration settings. This information can be generated upon
completion of a workflow with the command `snakemake -j1 summarize_methods`. This will create an additional output
file, `results/reports/methods_summary.md`, that contains the best description of the workflow as was actually run.
Note that this information focuses on methodology as opposed to results, and only requires the relevant conda
environments be installed; so if you want to predict what the workflow will do before actually running it,
complete user configuration, run `snakemake -j1 --use-conda --conda-create-envs-only`, and then run the `summarize_methods`
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
