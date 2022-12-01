library(testthat)

## include snakemake.s4 class for script interface
if (!isClass("snakemake.s4")) {
  source("snakemake.R")
}

source("aggregate_performance_metrics.R")

create.tmp.metrics.file <- function(parent.dirname, metrics.path) {
  dir.create(file.path(parent.dirname, dirname(metrics.path)), recursive = TRUE)
  out.fn <- file.path(parent.dirname, metrics.path)
  out.data <- data.frame(
    "s" = sample(1:100, 1),
    "h:m:s" = paste(sample(1:100, 1), sample(1:100, 1), sample(1:100, 1), sep = ":"),
    "max_rss" = sample(1:1000, 1),
    "max_vms" = sample(1:1000, 1),
    "max_ums" = sample(1:1000, 1),
    "max_pss" = sample(1:1000, 1),
    "io_in" = sample(1:1000, 1),
    "io_out" = sample(1:1000, 1),
    "mean_load" = sample(1:1000, 1),
    "cpu_time" = sample(1:1000, 1),
    check.names = FALSE
  )
  write.table(out.data, out.fn, row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
  out.data
}

test_that("aggregate.performance.metrics can load data from single metrics file", {
  in.dir <- tempdir()
  rulename <- "rulename1"
  expected <- create.tmp.metrics.file(in.dir, file.path(rulename, "subdir", "single_file1.tsv"))
  observed <- aggregate.performance.metrics(in.dir, rulename)
  expect_equal(observed, expected)
})

test_that("aggregate.performance.metrics can load data from multiple metrics files", {
  in.dir <- tempdir()
  rulename <- "rulename2"
  expected <- create.tmp.metrics.file(in.dir, file.path(rulename, "subdir1", "single_file1.tsv"))
  for (i in seq(2, 5)) {
    expected <- rbind(
      expected,
      create.tmp.metrics.file(in.dir, file.path(
        rulename,
        paste("subdir", i, sep = ""),
        paste("single_file", i, ".tsv", sep = "")
      ))
    )
  }
  observed <- aggregate.performance.metrics(in.dir, rulename)
  expect_equal(observed, expected)
})

test_that("aggregate.performance.metrics understands when a specified rule does not exist", {
  expect_error(aggregate.performance.metrics(tmpdir(), "rule1"))
})

test_that("aggregate.performance.metrics understands when the benchmark parent directory does not exist", {
  expect_error(aggregate.performance.metrics("FAKE_DIRNAME", "rule1"))
})
