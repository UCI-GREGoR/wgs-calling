library(testthat)

## include snakemake.s4 class for script interface
if (!isClass("snakemake.s4")) {
  source("snakemake.R")
}

source("construct_linker_from_inputs.R")

linker.df <- data.frame(
  subject = c(
    "internalid1", "internalid2",
    "internalid3", "internalid4",
    NA, "internalid6"
  ),
  jira = c(
    "RT-0001", "RT-0001",
    "RT-0002", "RT-0002",
    "RT-0003", NA
  ),
  ru = c(
    "RU00001", "RU00001",
    "RU00002", "RU00003",
    "RU00004", "RU00004"
  ),
  sq = c(
    "SQ0001", "SQ0002",
    NA, "SQ0004",
    "SQ0005", "SQ0006"
  ),
  ls = c(
    "LS00001", "LS00002",
    "LS00003", NA,
    "LS00005", "LS00006"
  ),
  sex = c(
    "Female", "Unknown",
    "Male", "Female",
    "Female", "NA"
  )
)

test_that("construct.output.stems can construct IDs from valid input data", {
  expected.out <- c(
    "internalid1_LS00001_SQ0001",
    "internalid2_LS00002_SQ0002",
    NA,
    NA,
    NA,
    "internalid6_LS00006_SQ0006"
  )
  observed.out <- construct.output.stems(linker.df)
  expect_identical(expected.out, observed.out)
})

test_that("construct.output.stems detects the absence of the required subject column", {
  broken.linker <- linker.df
  broken.linker$subject <- NULL
  expect_error(construct.output.stems(broken.linker))
})

test_that("construct.output.stems detects the absence of the required ls column", {
  broken.linker <- linker.df
  broken.linker$ls <- NULL
  expect_error(construct.output.stems(broken.linker))
})

test_that("construct.output.stems detects the absence of the required sq column", {
  broken.linker <- linker.df
  broken.linker$sq <- NULL
  expect_error(construct.output.stems(broken.linker))
})

test_that("apply.id.mappings recognizes ad hoc 'PMGRC' replacement patterns", {
  df <- data.frame(
    c(NA, "some filler text", "SAMPLE SWAP: this is actually PMGRC-123-456-7", "more filler text", NA),
    c("some filler text", "SAMPLE SWAP: this is actually PMGRC-987-654-3", NA, "more filler text", NA),
    c("some filler text", "more filler text", "even more filler text", "final filler text", NA),
    c(NA, NA, NA, NA, NA)
  )
  vec <- c("PMGRC-1-2-3", "PMGRC-4-5-6", "PMGRC-7-8-9", "PMGRC-11-22-3", "PMGRC-44-55-6")
  expected <- c(
    "PMGRC-1-2-3",
    "PMGRC-987-654-3",
    "PMGRC-123-456-7",
    "PMGRC-11-22-3",
    "PMGRC-44-55-6"
  )
  observed <- apply.id.mappings(df, vec)
  expect_identical(expected, observed)
})
