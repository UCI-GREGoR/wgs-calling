library(stringr, quietly = TRUE)

#' Combine subject ID, analyte ID, and sequencing index ID
#' into a '_'-delimited identifier for exported data preparation.
#'
#' @details
#' If any of the required columns have an NA entry for a particular subject,
#' the entire resultant ID is set to NA.
#'
#' @param df data.frame; input data minimally containing columns 'subject',
#' 'analyte', and 'index'
#' @return character vector; constructed output IDs for each row in input
#' data.frame, or NA if any of the required inputs for the row were NA
construct.output.stems <- function(df) {
  stopifnot(is.data.frame(df))
  stopifnot(
    "input data frame contains column 'subject'" = "subject" %in% colnames(df),
    "input data frame contains column 'analyte'" = "analyte" %in% colnames(df),
    "input data frame contains column 'index'" = "index" %in% colnames(df)
  )
  ## try to construct IDs of the format SUBJECTID_LSID_SQID
  output.stems <- paste(df$subject, df$analyte, df$index, sep = "_")
  output.stems[is.na(df$subject) | is.na(df$analyte) | is.na(df$index)] <- NA
  output.stems
}

#' Search for freetext ID remappings and apply them as needed
#'
#' The idea here is that there are evidently some freetext
#' notes along the lines of
#' "SAMPLE SWAP - this is actually PMGRC-\d+-\d+-\d" that
#' only exist in freetext in rows that are labeled by
#' _unmapped_ subject ID. For downstream analysis, these ID
#' mappings must be detected and applied. To make things
#' even more confusing, there doesn't appear to be any sense
#' to which columns of the logbook this information is reported in.
#'
#' The logic in this initial implementation will likely need to be
#' made more complex when different upstream freetext field types
#' are detected.
#'
#' @param df data.frame, input logbook sheet from read.xlsx
#' @param vec character vector, vector of actual subject IDs from
#' current sheet
#' @return character vector, mapped version of input character vector
#' with ID relabelings applied
apply.id.mappings <- function(df, vec) {
  stopifnot(is.data.frame(df))
  stopifnot(is.character(vec))
  stopifnot("input vector should represent subject ID column of data frame" = nrow(df) == length(vec))
  known.patterns <- c("^SAMPLE SWAP.*this is actually (PMGRC-\\d+-\\d+-\\d).*$" = "\\1")
  for (i in seq_len(ncol(df))) {
    new.ids <- rep(NA, length(vec))
    for (j in seq_len(length(known.patterns))) {
      has.pattern <- stringr::str_detect(df[, i], names(known.patterns)[j])
      has.pattern[is.na(has.pattern)] <- FALSE
      new.ids[has.pattern] <- stringr::str_replace(
        df[has.pattern, i],
        names(known.patterns)[j],
        known.patterns[j]
      )
    }
    vec[!is.na(new.ids)] <- new.ids[!is.na(new.ids)]
  }
  vec
}

#' Contextually overwrite or append linker information to an existing
#' data frame from a simple two column id->value linker file
#'
#' @details
#' If entirely new subjects are encountered in the linker file,
#' those subjects will be appended to the existing data frame
#' with relevant linker information in the appropriate column
#' and NA everywhere else
#'
#' @param df data.frame existing linker information, to which new
#' information will be overwritten or appended
#' @param linker.fn character; name of simple linker file with new
#' information to be added to the existing data frame
#' @param target.colname character; name of column in existing data
#' frame to which linker information will be saved. column does not
#' necessarily need to already exist, although it is intended to
#' @return data.frame; modified version of input data frame with
#' new information overwritten or appended
add.linker.data <- function(df, linker.fn, target.colname) {
  stopifnot(is.data.frame(df))
  stopifnot(
    is.character(linker.fn),
    length(linker.fn) == 1
  )
  stopifnot(is.character(target.colname))
  linker.df <- read.table(linker.fn, header = TRUE, sep = "\t", quote = "", comment.char = "", stringsAsFactors = FALSE)
  rownames(linker.df) <- linker.df[, 1]
  ## handle the situation where people not yet present are being added here
  new.subjects <- linker.df[!(linker.df[, 1] %in% df$subject), 1]
  if (length(new.subjects) > 0) {
    new.df <- data.frame(
      subject = new.subjects,
      jira = NA,
      project = NA,
      index = new.subjects,
      analyte = NA,
      sex = NA,
      external = NA
    )
    df <- rbind(df, new.df)
  }
  df[df[, "subject"] %in% linker.df[, 1], target.colname] <-
    linker.df[df[df[, "subject"] %in% linker.df[, 1], "subject"], 2]
  df
}

#' Run primary logic of this script, wrapped such that sourcing
#' this file will not cause actual code execution.
#'
#' @param sex.linker.fn character or NULL; name of input sex linker
#' file. Can be NULL, in which case it is effectively ignored
#' @param external.id.linker.fn character or NULL; name of input external
#' ID linker file. Can be NULL, in which case it is effectively ignored
#' @param out.fn character; name of file to which to write output linker data
run.construct.linker <- function(sex.linker.fn,
                                 external.id.linker.fn,
                                 out.fn) {
  stopifnot(
    is.character(sex.linker.fn) || is.null(sex.linker.fn),
    length(sex.linker.fn) <= 1
  )
  stopifnot(
    is.character(external.id.linker.fn) || is.null(external.id.linker.fn),
    length(external.id.linker.fn) <= 1
  )
  stopifnot(
    is.character(out.fn),
    length(out.fn) == 1
  )
  df <- data.frame(
    subject = c("A"),
    jira = c("A"),
    project = c("A"),
    index = c("A"),
    analyte = c("A"),
    sex = c("A")
  )
  df <- df[-1, ]
  if (!is.null(sex.linker.fn)) {
    df <- add.linker.data(df, sex.linker.fn, "sex")
  }
  if (!is.null(external.id.linker.fn)) {
    df <- add.linker.data(df, external.id.linker.fn, "external")
  }
  ## deal with possibility that no external ID has been provided
  df[is.na(df[, "external"]), "external"] <- df[is.na(df[, "external"]), "subject"]

  ## deal with apparent internal newlines
  for (index in seq_len(ncol(df))) {
    df[, index] <- stringr::str_replace_all(df[, index], "\n", " ")
  }

  ## make encodings of missing sex be "Unknown" instead of possibly NA
  df[is.na(df[, "sex"]), "sex"] <- "Unknown"

  write.table(df, out.fn, row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
}

if (exists("snakemake")) {
  run.construct.linker(
    snakemake@params[["sex_linker"]],
    snakemake@params[["external_id_linker"]],
    snakemake@output[["linker"]]
  )
}
