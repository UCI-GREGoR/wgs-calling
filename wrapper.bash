#!/usr/bin/env bash
snakemake -j50 --profile ../sge-profile -p --rerun-incomplete --use-conda --rerun-triggers mtime
