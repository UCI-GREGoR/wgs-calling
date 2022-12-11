#' Create a placeholder snakemake S4 class
#' for high-fidelity unit testing of
#' embedded R scripts

#' Function to verify consistency of object
#' of S4 snakemake class. The only meaningful
#' restriction beyond the class definition
#' is that the list entries should be
#' character vectors.
#'
#' @param object snakemake.s4 object to test
#' @return if the object is valid, the
#' logical TRUE; else a character vector
#' of reasons the object is non-compliant.
check.snakemake.s4 <- function(object) {
  errors <- c()
  if (!all(sapply(object@input, is.character))) {
    errors <- c(errors, "input list contains non-character entries")
  }
  if (!all(sapply(object@output, is.character))) {
    errors <- c(errors, "output list contains non-character entries")
  }
  if (length(object@output) == 0) {
    errors <- c(errors, "output list needs to contain at least one target")
  }
  ifelse(length(errors) == 0, TRUE, errors)
}

#' snakemake S4 class
#'
#' note that params has some odd conversion behaviors from
#' whatever is specified in python.
#'
#' - None in the params block becomes NULL in R
#'
#' @param input list of input files; should all be character,
#' and can be empty
#' @param output list of output files; should be all character,
#' but must have at least length 1
#' @param params list of params entries; no specific type restrictions
#'
setClass("snakemake.s4", representation(
  input = "list",
  output = "list",
  params = "list"
), validity = check.snakemake.s4)
