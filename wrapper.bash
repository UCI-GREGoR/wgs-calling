#!/usr/bin/env bash
snakemake -j50 --profile ../sge-profile -p --rerun-incomplete --rerun-triggers mtime --use-conda
