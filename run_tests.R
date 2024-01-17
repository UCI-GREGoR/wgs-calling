#!/usr/bin/env Rscript
library(testthat)
library(covr)
setwd("workflow/scripts")
target.files <- list.files(".", "^[^\\-]*R$")
target.files <- target.files[target.files != "snakemake.R"]
test.files <- list.files(".", "^test-.*R$")
res <- covr::file_coverage(target.files, test.files)
covr:::to_codecov(res)
