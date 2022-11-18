#!/usr/bin/env bash
snakemake -j50 --profile ../sge-profile -p --rerun-incomplete --use-conda -F --rerun-triggers mtime
