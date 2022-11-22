# Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- legacy support for fastQC, fastp.
- legacy support for bwa-mem2 for alignment.
- legacy support for markdups.
- legacy support for aligned read QC (mostly picard).

### Changed

- `bwa-mem2` is pulled from conda, which currently pulls v2.2.1. the docker image
  from the legacy workflow has pulled from an untagged feature branch. this is bad.
- markdups switched from `samtools markdup` to `gatk MarkDuplicates` (i.e. picard).

### Deprecated

- `fastp` emits output fastqs with adapters trimmed, but these trimmed fastqs are not passed to `bwa-mem2`;
  this is for concordance with legacy behavior. this behavior is flagged for removal; adapter trimming is good.

[//]: # (- Added)
[//]: # (- Changed)
[//]: # (- Deprecated)
[//]: # (- Removed)
[//]: # (- Fixed)
[//]: # (- Security)
