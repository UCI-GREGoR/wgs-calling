#' For a rule, gather all available performance
#' benchmark data in the workflow installation.
#'
#' @param parent.dir character vector; path to benchmark
#' data parent directory. this will be something like
#' `results/performance_benchmarks`.
#' @param rule character vector; name of rule for which
#' to aggregate performance metrics. the rule name should
#' be a depth 1 subdirectory of the parent benchmark directory.
#' @return data.frame; aggregated observations from
#' each instance of the target rule in tabular format.
aggregate.performance.metrics <- function(parent.dir, rule) {
  ## input consistency checks
  stopifnot(is.character(parent.dir))
  stopifnot(length(parent.dir) == 1)
  stopifnot(is.character(rule))
  stopifnot(length(rule) == 1)
  stopifnot(dir.exists(parent.dir))
  rule.dir <- file.path(parent.dir, rule)
  stopifnot(dir.exists(rule.dir))

  ## aggregate all available performance metrics for a rule.
  ## note that artifacts from prior runs will persist, and so
  ## it is possible to end up with situations where metrics
  ## from old runs are included in new runs' summaries.
  benchmark.files <- list.files(rule.dir, "*.tsv", recursive = TRUE, full.names = TRUE)

  ## the output is snakemake's aggregated psutils output. per
  ## https://stackoverflow.com/questions/46813371/meaning-of-the-benchmark-variables-in-snakemake,
  ## they will all have the format:
  ## s - running time in seconds
  ## h:m:s - running time in hour:minute:second format
  ## max_rss - maximum resident set size, non-swapped physical memory used by process
  ## max_vms - maximum virtual memory size, total virtual memory used by process
  ## max_uss - maximum unique set size, memory unique to a process which would be freed by killing job
  ## max_pss - maximum proportional set size
  ## io_in - cumulative amount of MB read in
  ## io_out - cumulative amount of MB written out
  ## mean_load - CPU usage over time, divided by total runtime (in seconds)
  ## cpu_time - CPU time summed for user and system
  res <- data.frame()
  for (benchmark.file in benchmark.files) {
    df <- read.table(benchmark.file, comment.char = "", checkNames = FALSE, quote = "", stringsAsFactors = FALSE)
    if (nrow(res) == 0) {
      res <- df
    } else {
      res <- rbind(res, df)
    }
  }
  res
}
