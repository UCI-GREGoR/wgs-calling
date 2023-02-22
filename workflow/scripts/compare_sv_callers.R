#' Store preset ggplot theme data for use with all relevant output plots
my.theme <- theme_light() + theme(
  plot.title = element_text(size = 16, hjust = 0.5),
  axis.title = element_text(size = 14),
  axis.text = element_text(size = 12),
  strip.background = element_blank(),
  strip.text = element_text(size = 14, colour = "black")
)

#' Load sv source information from file and combine across sample.
#' Annotate data frame with source information per caller, with the
#' intention of using this information for upset.
#'
#' @param filenames character vector; input filenames containing summary
#' data from bcftools query
#' @param sv.callers character vector; names of requested sv calling tools
#' provided by user in configuration
#' @return data.frame; combined input data with harmonized headers and columns
#' indicating whether the variant (or something very similar) was
#' found in the output of each of the sv callers, as well as how many variants from that
#' caller contributed to the final variant
load.data <- function(filenames, sv.callers) {
  res <- data.frame()
  for (filename in filenames) {
    h <- read.table(filename,
      header = FALSE, stringsAsFactors = FALSE, sep = "\t",
      comment.char = "", quote = ""
    )
    if (nrow(res) == 0) {
      res <- h
    } else {
      res <- rbind(res, h)
    }
  }
  if (length(colnames(res)) == 2) {
    colnames(res) <- c("SVTYPE", "variant_origin")
  } else {
    colnames(res) <- c(
      "CHROM", "POS", "ID", "REF", "ALT",
      "QUAL", "FILTER", "SVTYPE", "variant_origin"
    )
  }
  ## So, there are two questions one can ask about source information for each output variant.
  ## There's the question of whether a particular caller contributed to a variant at all;
  ## and then in addition there's the question of how many variants were collapsed into a single
  ## variant from each caller output file. svdb seems to report a non-unique list of source files
  ## in INFO/svdb_origin, and so both questions can be addressed, if messily.
  ## 
  ## Since writing the above, truvari support has been added, and this tracking information is
  ## coming from a somewhat different source. The column header has been updated to reflect this
  ## more generic origin.
  for (sv.caller in sv.callers) {
    res[, paste("in", sv.caller, sep = ".")] <- grepl(sv.caller, res[, "variant_origin"])
    res[, paste("count", sv.caller, sep = ".")] <- stringr::str_count(res[, "variant_origin"], sv.caller)
  }
  res
}

#' Using an input data.frame, generate a list of source
#' information compatible with upset
#'
#' @param df data.frame; loaded input data from load.data
#' @param sv.callers character vector; names of requested sv calling tools
#' provided by user in configuration
#' @return list; sv caller source information compatible with
#' upset fromList
generate.upset.list <- function(df, sv.callers) {
  res <- list()
  for (sv.caller in sv.callers) {
    res[[sv.caller]] <- which(df[, paste("in", sv.caller, sep = ".")])
  }
  res
}

#' Based on a target set of SV callers, generate counts of variants
#' per SV type, for use in later barplots
#'
#' @param df data.frame; loaded input data from load.data
#' @param sv.callers character vector; names of requested sv calling tools
#' provided by user in configuration
#' @return data.frame; tidy barplot input for all sv callers
generate.barplot.input <- function(df, sv.callers) {
  caller <- c()
  sv <- c()
  sv.count <- c()
  for (sv.caller in sv.callers) {
    for (var.type in unique(df$SVTYPE)) {
      caller <- c(caller, sv.caller)
      sv <- c(sv, var.type)
      sv.count <- c(sv.count, sum(df[df$SVTYPE == var.type, paste("count", sv.caller, sep = ".")]))
    }
  }
  res <- data.frame(
    sv.caller = caller,
    sv.type = sv,
    sv.count = sv.count
  )
  res
}

#' Using an input data.frame, create a tool-specific barplot
#' representing the number of variants of each SVTYPE that were
#' present in the initial callset
#'
#' @param df data.frame; prepared barplot input data from generate.barplot.input
#' @param sv.caller character vector; name of requested sv calling tool
#' @return ggplot2 barplot
get.tool.variant.barplot <- function(df, sv.caller) {
  ## compute tidy format data frame for barplot, based on variant
  ## types observed by SV caller
  plot.data <- df[df$sv.caller == sv.caller, ]
  my.plot <- ggplot(aes(x = sv.type, y = sv.count), data = plot.data) + my.theme
  my.plot <- my.plot + geom_bar(stat = "identity")
  my.plot <- my.plot + xlab("SV Type") + ylab("Number of Variants of Type in Initial Calls")
  my.plot <- my.plot + scale_y_continuous(limits = c(0, max(df$sv.count * 1.2)))
  my.plot
}
