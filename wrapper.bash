#!/usr/bin/env bash

## stability update: there are a ton of conda environments in the full run of this workflow.
## during the initial run of a workflow, the DAG builds, remotes are connected, and then
## conda environments are linearly constructed, very slowly. this creates situations where
## the remotes then fail after the environments are constructed. at the same time, constructing
## the DAG itself can be a costly operation for some runs.
## as a hack: if no conda environments are present at all, run a preflight pass that
## exclusively constructs conda environments.
if [[ ! -d .snakemake/conda ]] ; then
    snakemake -j1 -p --rerun-incomplete --rerun-triggers mtime --use-conda --use-singularity --conda-create-envs-only
fi
snakemake -j150 --profile ../sge-profile -p --rerun-incomplete --rerun-triggers mtime --use-conda --use-singularity --cluster-config config/cluster.yaml
