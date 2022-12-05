library(testthat)

## include snakemake.s4 class for script interface
if (!isClass("snakemake.s4")) {
  source("snakemake.R")
}

source("aggregate_performance_metrics.R")

demo.data <- data.frame(
  "s" = sample(1:100, 1),
  "h:m:s" = paste(sample(1:100, 1), sample(1:100, 1), sample(1:100, 1), sep = ":"),
  "max_rss" = sample(1:1000, 1),
  "max_vms" = sample(1:1000, 1),
  "max_uss" = sample(1:1000, 1),
  "max_pss" = sample(1:1000, 1),
  "io_in" = sample(1:1000, 1),
  "io_out" = sample(1:1000, 1),
  "mean_load" = sample(1:1000, 1),
  "cpu_time" = sample(1:1000, 1),
  check.names = FALSE
)

create.tmp.metrics.file <- function(parent.dirname, metrics.path) {
  dir.create(file.path(parent.dirname, dirname(metrics.path)), recursive = TRUE)
  out.fn <- file.path(parent.dirname, metrics.path)
  out.data <- demo.data
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
  expect_null(aggregate.performance.metrics(tempdir(), "rule1"))
})

test_that("aggregate.performance.metrics understands when the benchmark parent directory does not exist", {
  expect_error(aggregate.performance.metrics("FAKE_DIRNAME", "rule1"))
})

test_that("get.io.description returns a seemingly-relevant character vector", {
  output <- get.io.description()
  expect_true(stringr::str_detect(output, stringr::regex("^\n\nInput and output represent.*\n\n$", dotall = TRUE)))
})

test_that("get.memory.description returns a seemingly-relevant character vector", {
  output <- get.memory.description()
  expect_true(stringr::str_detect(output, stringr::regex(paste("^\n\nRSS, resident set size.*",
    "VMS, virtual memory size.*",
    "USS, unique set size.*",
    "PSS, proportional set size.*\n\n$",
    sep = ""
  ), dotall = TRUE)))
})

test_that("get.cpu.walltime.description returns a seemingly-relevant character vector", {
  output <- get.cpu.walltime.description()
  expect_true(stringr::str_detect(output, stringr::regex("^\n\nCPU/walltime is a representation.*\n\n$",
    dotall = TRUE
  )))
})

test_that("get.walltime.description returns a seemingly-relevant character vector", {
  output <- get.walltime.description()
  expect_true(stringr::str_detect(output, stringr::regex("^\n\nWalltime represents.*\n\n$", dotall = TRUE)))
})

## testing ggplot functions is an imperfect art
test_that("walltime.plot returns a minimally-conformant plot", {
  output <- walltime.plot(demo.data, "rulename")
  ## check that standard theme has been applied
  expect_equal(output$theme$plot.title$hjust, 0.5)
  ## check that it is geom_point
  expect_true(!is.null(output$layers))
  expect_true(inherits(output$layers[[1]]$geom, "GeomPoint"))
  ## check that it just has the one layer
  expect_equal(length(output$layers), 1)
})

test_that("cpu.div.walltime.plot returns a minimally-conformant plot", {
  output <- cpu.div.walltime.plot(demo.data, "rulename", 2)
  ## check that standard theme has been applied
  expect_equal(output$theme$plot.title$hjust, 0.5)
  ## check that it is geom_point
  expect_true(!is.null(output$layers))
  expect_true(inherits(output$layers[[1]]$geom, "GeomPoint"))
  ## check that it has two layers: the points and the threads line
  expect_equal(length(output$layers), 2)
  expect_true(inherits(output$layers[[2]]$geom, "GeomHline"))
})

test_that("memory.plot returns a minimally-conformant plot", {
  output <- memory.plot(demo.data)
  ## check that standard theme has been applied
  expect_equal(output$theme$plot.title$hjust, 0.5)
  ## check that it is geom_point
  expect_true(!is.null(output$layers))
  expect_true(inherits(output$layers[[1]]$geom, "GeomPoint"))
  ## check that it just has the one layer
  expect_equal(length(output$layers), 1)
})

test_that("io.plot returns a minimally-conformant plot", {
  output <- io.plot(demo.data)
  ## check that standard theme has been applied
  expect_equal(output$theme$plot.title$hjust, 0.5)
  ## check that it is geom_point
  expect_true(!is.null(output$layers))
  expect_true(inherits(output$layers[[1]]$geom, "GeomPoint"))
  ## check that it just has the one layer
  expect_equal(length(output$layers), 1)
})
