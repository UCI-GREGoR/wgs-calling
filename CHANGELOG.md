# Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased] - DeepVariant

First support for DeepVariant. High burden tests are still running, and the UGE HPC profile
is very recently changed to more strictly enforce resource requests, so everything is in flux.

### Added

- support for DeepVariant
  - specify with `config/config.yaml` `snv-caller: deepvariant`
- automatically record performance benchmark data
  - stored under `results/performance_benchmarks`
- optionally create performance benchmark report for review during optimization
  - the report is not created by default, and requires substantially more work
    before it's really ready for primetime
  - to create the report, run `snakemake -j1 --use-conda performance_benchmarks`

## [0.1.0] - 2022-11-28

Integration of initial implementation branch, so that default has things of interest for casual viewing.
The initial high-burden test of the implementation is still running, so rapid and immediate changes
are anticipated.

### Added

- legacy support for fastQC, fastp.
- legacy support for bwa-mem2 for alignment.
- legacy support for markdups.
- legacy support for aligned read QC (mostly picard).
- legacy support for octopus.
- legacy support for tiddit, manta SV calling.
- legacy support for duphold.
- legacy support for merging tiddit/manta output.
- test coverage for all snakemake lambda functions (under `lib/`).
- test coverage for all added internal python scripts (under `workflow/scripts/`).

### Changed

- `bwa-mem2` is pulled from conda, which currently pulls v2.2.1. the docker image
  from the legacy workflow has pulled from an untagged feature branch. this is bad.
- markdups switched from `samtools markdup` to `gatk MarkDuplicates` (i.e. picard).
- `octopus` is orphaned, and the available versions are of questionable desirability.




### Deprecated

- `fastp` emits output fastqs with adapters trimmed, but these trimmed fastqs were not passed to `bwa-mem2` in legacy.
  this behavior should no longer be the default, but can be turned back on with the corresponding config/behaviors setting.

### Removed

- some filters on the tiddit output in particular were undocumented and are at the least temporarily removed,
  until such time as I can confirm what they were intended to do.

### Fixed

- `duphold` is actually being applied now. previously, it was only, and somewhat wonkily,
  being applied to tiddit.


[//]: # [Unreleased]

[//]: # (- Added)
[//]: # (- Changed)
[//]: # (- Deprecated)
[//]: # (- Removed)
[//]: # (- Fixed)
[//]: # (- Security)
