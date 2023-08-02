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
    "JT0001", "JT0001",
    "JT0002", "JT0002",
    "JT0003", NA
  ),
  project = c(
    "PR00001", "PR00001",
    "PR00002", "PR00003",
    "PR00004", "PR00004"
  ),
  index = c(
    "SI0001", "SI0002",
    NA, "SI0004",
    "SI0005", "SI0006"
  ),
  analyte = c(
    "AN00001", "AN00002",
    "AN00003", NA,
    "AN00005", "AN00006"
  ),
  sex = c(
    "Female", "Unknown",
    "Male", "Female",
    "Female", "NA"
  )
)

test_that("construct.output.stems can construct IDs from valid input data", {
  expected.out <- c(
    "internalid1_AN00001_SI0001",
    "internalid2_AN00002_SI0002",
    NA,
    NA,
    NA,
    "internalid6_AN00006_SI0006"
  )
  observed.out <- construct.output.stems(linker.df)
  expect_identical(expected.out, observed.out)
})

test_that("construct.output.stems detects the absence of the required subject column", {
  broken.linker <- linker.df
  broken.linker$subject <- NULL
  expect_error(construct.output.stems(broken.linker))
})

test_that("construct.output.stems detects the absence of the required analyte column", {
  broken.linker <- linker.df
  broken.linker$analyte <- NULL
  expect_error(construct.output.stems(broken.linker))
})

test_that("construct.output.stems detects the absence of the required index column", {
  broken.linker <- linker.df
  broken.linker$index <- NULL
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

test_that("add.linker.data correctly replaces and appends data as appropriate", {

})
