# Quick Overview for Pipeline Configuration

Note that this is a brief runthrough of rapid configuration
of the workflow for deployment. Additional information for
finer-grained configuration is available in [the readme](README.md).

## One time only: install mamba
- Follow instructions [here](https://mamba.readthedocs.io/en/latest/installation.html) to install mambaforge.
- Create a deployment environment containing snakemake:
    - `mamba install -c bioconda -c conda-forge -n snakemake "snakemake>=7"`
- Download a suitable HPC scheduler profile.
    - slurm support is now (snakemake 7.25 or so) cooked into snakemake.
    - profiles are available [here](https://github.com/Snakemake-Profiles).
    - for internal UGE support, a pre-configured profile [is available](git@gitlab.com:lightning.auriga1/pmgrc-sge-profile.git).

## Configuration for each workflow deployment
- Clone repo from [this github/ssh link](git@gitlab.com:lightning.auriga1/wgs-pipeline.git).
- Checkout latest release version.
- Update configuration for run in `config/config.yaml`:
    - choose a run mode in `behaviors/outcome`. Options are `fastq`, `alignment`, `calling`, and `release`.
      For most runs, setting this first to `calling` is appropriate, as it will hopefully flag any issues
      with the data and provide a point where the user can inspect QC results before proceeding. Once the
      run has completed to the user's satisfaction, set this to `release` and rerun, and it will prepare
      a data export (cram and crai for alignments; vcf.gz and tbi for SNV and SV calls; and checksums)
      and send it to one of the configured export locations (see below). Alternatively, you can set this
      to `release` from the start and live on the edge.
    - for most analyses, you will need to adjust the following, but only if their functionality is required:
        - `behaviors/import-s3/profile-name`: name of aws-cli profile in environment, for accessing remote s3 input fastqs.
          If you have local fastqs, this is not required. Remote fastqs are detected as having `s3://` prefixes
          in the input manifest.
        - `behaviors/export-directory`: path to local directory for data export. This is only required if you're running
          in `release` mode and want a local data export copy.
        - `behaviors/export-s3/bucket-name` and `behaviors/export-s3/profile-name`: name of s3 bucket to which to send copy
          of final data, and name of aws-cli profile in the environment to use for authentication. As above,
          only required if running in `release` mode and want a remote data export copy.
- Update manifest for run in `config/manifest.tsv`:
    - from fastq input:
        - each row corresponds to a fastq pair.
        - the default manifest in the repo has the correct column headers in place.
        - fastqs R1 and R2 should be provided as absolute or relative paths, or `s3://` links. The prefix
          `s3://` is detected and triggers a download from aws. If required, set the configuration setting
          `behaviors/import-s3/profile-name`, described above.
        - lanes should be `L001`, `L002`, etc.; or `combined` if only a single fastq per-sample is provided.
        - `projectid` is for legacy support and likely to be deprecated; it was meant to correspond to flowcell,
          such that the workflow would support running multiple flowcells in the same run but split out by meaningful
          batch. Until this functionality is removed, set this to some constant alphanumeric value.
    - from bam (pre-aligned read) input:
        - each row corresponds to a single bam.
        - the default manifest should be modified to have headers `projectid`, `sampleid`, and `bam`.
        - see above for description of `projectid`.
        - bams should be provided as absolute or relative paths. Remotes not currently supported.
- Provide sample descriptors:
    - `config/sample_linking_external_ids.tsv`:
        - first column is `sampleid`, corresponding to entry in manifest.
        - second column is `externalid`, to which the sample ID will be mapped in prepared output files.
    - `config/sample_linking_sex.tsv`:
        - first column is `sampleid`, corresponding to entry in manifest.
        - second column is `sex`, sample self-reported sex. Values should be `Female` or `Male`, case sensitive.
          Samples without linker data in this file will be set to missing. This information is used exclusively
          for sexcheck in [somalier](https://github.com/brentp/somalier).
- Update `wrapper.bash` to point to the installed UGE profile directory.

## Running the workflow

- Launch the workflow: `qsub -N run_pipeline -V -q small -cwd -S /bin/bash ./wrapper.bash`, for example.
- Wait.
- Check UGE job status.
    - `qstat` shows a tabular summary of running jobs, for example.
    - workflow run length varies wildly depending on input size and configuration settings.
    - expect at least 1-2 days for full run from scratch.
- Check workflow output, in `run_pipeline.e*` if the above command was used.
    - a nice simple check is `grep Error run_pipeline.e*` (note case). If it ran successfully, this should return nothing.
    - `tail run_pipeline.e*` should contain lines like:
        - `3000 of 3000 steps (100%) done`
        - `Complete log: .snakemake/log/{some stuff}.snakemake.log`
- If it goes down with `Error`s:
    - debug.
    - file bug report, if developer is available.
    - just try rerunning as-is, sometimes HPCs are unstable.
- Check multiqc output, at `results/multiqc/{project_id}/*.html`.
- Good luck!
