#' For convenience, define a single ggplot theme/style
my.theme <- theme_light() + theme(
  plot.title = element_text(size = 16, hjust = 0.5),
  axis.title = element_text(size = 14),
  axis.text = element_text(size = 12),
  strip.background = element_blank(),
  strip.text = element_text(size = 14, colour = "black"),
  legend.title = element_text(size = 14),
  legend.text = element_text(size = 13)
)

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
    df <- read.table(benchmark.file,
      comment.char = "", check.names = FALSE, quote = "",
      stringsAsFactors = FALSE, header = TRUE
    )
    if (nrow(res) == 0) {
      res <- df
    } else {
      res <- rbind(res, df)
    }
  }
  res
}

#' Create plot of distribution of simple walltime in hours
#'
#' @param benchmark.df data.frame; input aggregated rule
#' benchmarking plot from snakemake benchmark directive
#' @param rule.name character vector; name of relevant
#' rule corresponding to these metrics
#' @return ggplot object
walltime.plot <- function(benchmark.df, rule.name) {
  benchmark.df$rule.name <- rule.name
  benchmark.df$s <- benchmark.df$s / 3600
  max.y <- max(benchmark.df$s) * 1.2
  my.plot <- ggplot(aes(x = rule.name, y = s), data = benchmark.df)
  my.plot <- my.plot + my.theme + geom_point(position = "jitter")
  my.plot <- my.plot + xlab("") + ylab("Walltime, in hours")
  my.plot <- my.plot + scale_y_continuous(limits = c(0, max.y))
  my.plot
}

#' Create plot of distribution of CPU/walltime ratio
#'
#' @param benchmark.df data.frame; input aggregated rule
#' benchmarking plot from snakemake benchmark directive
#' @param rule.name character vector; name of relevant
#' rule corresponding to these metrics
#' @param provided.threads integer; number of threads
#' configured for the particular rule, which should technically
#' be the theoretical maximum of the computed ratio
#' @return ggplot object
cpu.div.walltime.plot <- function(benchmark.df, rule.name, provided.threads) {
  benchmark.df$rule.name <- rule.name
  benchmark.df$cpu.walltime.ratio <- benchmark.df[, "cpu_time"] / benchmark.df$s
  max.y <- max(c(benchmark.df$cpu.walltime.ratio, provided.threads)) * 1.2
  my.plot <- ggplot(aes(x = rule.name, y = cpu.walltime.ratio), data = benchmark.df)
  my.plot <- my.plot + my.theme + geom_point(position = "jitter")
  my.plot <- my.plot + xlab("") + ylab("CPU/Walltime")
  my.plot <- my.plot + scale_y_continuous(limits = c(0, max.y))
  if (!is.na(provided.threads)) {
    my.plot <- my.plot + geom_hline(yintercept = provided.threads, lty = 2)
  }
  my.plot
}

#' Create plot of distribution of memory allocation in MB
#'
#' @param benchmark.df data.frame; input aggregated rule
#' benchmarking plot from snakemake benchmark directive
#' @return ggplot object
memory.plot <- function(benchmark.df) {
  mem.types <- c(
    "max_rss" = "Max RSS",
    "max_vms" = "Max VMS",
    "max_uss" = "Max USS",
    "max_pss" = "Max PSS"
  )
  plot.data <- data.frame(
    x = rep(mem.types, each = nrow(benchmark.df)),
    y = c(
      benchmark.df[, names(mem.types)[1]],
      benchmark.df[, names(mem.types)[2]],
      benchmark.df[, names(mem.types)[3]],
      benchmark.df[, names(mem.types)[4]]
    )
  )
  plot.data$x <- factor(plot.data$x, levels = mem.types)
  max.y <- max(plot.data$y) * 1.2
  my.plot <- ggplot(aes(x = x, y = y), data = plot.data)
  my.plot <- my.plot + my.theme + geom_point(position = "jitter")
  my.plot <- my.plot + xlab("") + ylab("Memory Used, in MB")
  my.plot <- my.plot + scale_y_continuous(limits = c(0, max.y))
  my.plot
}

#' Create plot of I/O totals for tasks
#'
#' @param benchmark.df data.frame; input aggregated rule
#' benchmarking plot from snakemake benchmark directive
#' @return ggplot object
io.plot <- function(benchmark.df) {
  io.types <- c(
    "io_in" = "Input",
    "io_out" = "Output"
  )
  plot.data <- data.frame(
    x = rep(io.types, each = nrow(benchmark.df)),
    y = c(
      benchmark.df[, names(io.types)[1]],
      benchmark.df[, names(io.types)[2]]
    )
  )
  plot.data$x <- factor(plot.data$x, levels = io.types)
  max.y <- max(plot.data$y) * 1.2
  my.plot <- ggplot(aes(x = x, y = y), data = plot.data)
  my.plot <- my.plot + my.theme + geom_point(position = "jitter")
  my.plot <- my.plot + xlab("") + ylab("Data Transferred, in MB")
  my.plot <- my.plot + scale_y_continuous(limits = c(0, max.y))
  my.plot
}

#' Get a human-readable description of the benchmarking
#' walltime metric
#'
#' @return character vector; intended for emission in markdown
#' with cat()
get.walltime.description <- function() {
  "\n\nWalltime represents the actual amount of time consumed by the entire job, end to end.\n\n"
}

#' Get a human-readable description of the benchmarking
#' cpu/walltime metric
#'
#' @return character vector; intended for emission in markdown
#' with cat()
get.cpu.walltime.description <- function() {
  paste("\n\nCPU/walltime is a representation of the amount of possible CPU activity that is ",
    "actually achieved by a process. If a job is allocated N threads, a ratio of CPU/walltime = N ",
    "means that all threads were saturated with tasks during runtime. Lower numbers mean that ",
    "some bottleneck other than direct CPU processivity is blocking your job; often that is I/O speed, ",
    "but it can also indicate large imbalances in how tasks are apportioned to threads.\n\n",
    sep = ""
  )
}

#' Get a human-readable description of the benchmarking
#' memory metrics
#'
#' @return character vector; intended for emission in markdown
#' with cat()
get.memory.description <- function() {
  paste("\n\nRSS, resident set size, is the amount of non-swapped physical memory a process has used. ",
    "VMS, virtual memory size, is the total amount of virtual memory used by a process. ",
    "USS, unique set size, is the amount of memory unique to a particular process. ",
    "PSS, proportional set size, is the amount of memory shared with other processes, balanced across ",
    "all processes. All metrics are presented as max values across the total process runtime.\n\n",
    sep = ""
  )
}

#' Get a human-readable description of the benchmarking
#' I/O metrics
#'
#' @return character vector; intended for emission in markdown
#' with cat()
get.io.description <- function() {
  paste("\n\nInput and output represents the amount of data transferred from disk or network ",
    "by a process. Large I/O values combined with a CPU/walltime ratio less than N(threads) ",
    "may suggest potential for optimization by transferring a task to a system with ",
    "network optimization.\n\n",
    sep = ""
  )
}
