#' Load sv source information from file and combine across sample.
#' Annotate data frame with source information per caller, with the
#' intention of using this information for upset.
#'
#' @param filenames character vector; input filenames containing summary
#' data from bcftools query
#' @param sv.callers character vector; names of requested sv calling tools
#' provided by user in configuration
#' @return data.frame; combined input data with harmonized headers and logical
#' vector columns indicating whether the variant (or something very similar) was
#' found in the output of each of the sv callers
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
  colnames(res) <- c(
    "CHROM", "POS", "ID", "REF", "ALT",
    "QUAL", "FILTER", "SVTYPE", "svdb_origin"
  )
  for (sv.caller in sv.callers) {
    res[, paste("in", sv.caller, sep = ".")] <- grepl(sv.caller, res[, "svdb_origin"])
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
