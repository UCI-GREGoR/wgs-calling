library(openxlsx, quietly = TRUE)
library(stringr, quietly = TRUE)

#' Combine subject ID, ls (analyte) ID, and sq (sequencing index) ID
#' into a '_'-delimited identifier for exported data preparation.
#'
#' @details
#' If any of the required columns have an NA entry for a particular subject,
#' the entire resultant ID is set to NA.
#'
#' @param df data.frame; input data minimally containing columns 'subject',
#' 'ls', and 'sq'
#' @return character vector; constructed output IDs for each row in input
#' data.frame, or NA if any of the required inputs for the row were NA
construct.output.stems <- function(df) {
  stopifnot(is.data.frame(df))
  stopifnot(
    "input data frame contains column 'subject'" = "subject" %in% colnames(df),
    "input data frame contains column 'ls'" = "ls" %in% colnames(df),
    "input data frame contains column 'sq'" = "sq" %in% colnames(df)
  )
  ## try to construct IDs of the format SUBJECTID_LSID_SQID
  output.stems <- paste(df$subject, df$ls, df$sq, sep = "_")
  output.stems[is.na(df$subject) | is.na(df$ls) | is.na(df$sq)] <- NA
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

parse.logbook <- function(input.fn, output.fn) {
  stopifnot(
    is.character(input.fn),
    length(input.fn) == 1
  )
  stopifnot(
    is.character(output.fn),
    length(output.fn) == 1
  )
  sheet.names <- openxlsx::getSheetNames(input.fn)
  subject.id <- c()
  jira.tickets <- c()
  sq.id <- c()
  ls.id <- c()
  ru.id <- c()
  sex.id <- c()
  for (sheet.name in sheet.names) {
    df <- openxlsx::read.xlsx(input.fn, sheet = sheet.name, check.names = FALSE)
    colnames(df) <- tolower(colnames(df))
    if (colnames(df)[1] == "subject.id") {
      ## resolve chaos
      subject.col <- 1
      jira.col <- which(stringr::str_detect(colnames(df), "jira\\.ticket\\.for\\.batches\\.in\\.flight"))
      sq.col <- which(stringr::str_detect(colnames(df), "sq\\.id"))
      ls.col <- which(stringr::str_detect(colnames(df), "ls\\.id"))
      ru.col <- which(stringr::str_detect(colnames(df), "ru\\.id"))
      sex.col <- which(stringr::str_detect(colnames(df), "biological\\.sex"))
      if (length(jira.col) == 0) {
        jira.col <- which(stringr::str_detect(colnames(df), "^jira\\.ticket\\(s\\)\\.\\("))
      }
      if (length(jira.col) != 1) {
        next
      }
      ## "fix": upstream logbook has sporadic freetext annotations that indicate when a
      ## subject needs to have an ID update applied to it. for some IDs, this information
      ## is not recorded but rather just applied without annotation. It's not clear what
      ## governs the difference between those two situations
      df[, subject.col] <- apply.id.mappings(df, df[, subject.col])
      for (row.num in seq_len(nrow(df))) {
        subject.id <- c(subject.id, df[row.num, subject.col])
        jira.tickets <- c(jira.tickets, df[row.num, jira.col])
        sq.id <- c(sq.id, df[row.num, sq.col])
        ls.id <- c(ls.id, df[row.num, ls.col])
        ru.id <- c(ru.id, df[row.num, ru.col])
        if (length(sex.col) == 1) {
          sex.id <- c(sex.id, df[row.num, sex.col])
        } else {
          sex.id <- c(sex.id, "Unknown")
        }
      }
    }
  }
  res <- data.frame(
    subject = subject.id,
    jira = jira.tickets,
    ru = ru.id,
    sq = sq.id,
    ls = ls.id,
    sex = sex.id
  )
  ## remove entries with useless content
  res <- res[!(is.na(res[, 1]) & is.na(res[, 2]) & is.na(res[, 3]) & is.na(res[, 4]) & is.na(res[, 5])), ]
  ## construct output stems for deliverables, if sufficient information is present to fit the accepted format
  output.stems <- construct.output.stems(res)
  res$output <- output.stems
  res
}

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
      ru = NA,
      sq = NA,
      ls = NA,
      sex = NA,
      output = NA
    )
    df <- rbind(df, new.df)
  }
  df[, target.colname] <- linker.df[df[, "subject"], 2]
  df
}

run.construct.linker <- function(logbook.fn,
                                 sex.linker.fn,
                                 external.id.linker.fn,
                                 out.fn) {
  stopifnot(
    is.character(logbook.fn) || is.null(logbook.fn),
    length(logbook.fn) <= 1
  )
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
    ru <- c("A"),
    sq <- c("A"),
    ls <- c("A"),
    sex <- c("A")
  )
  df <- df[-1, ]
  if (!is.null(logbook.fn)) {
    df <- parse.logbook(logbook.fn, out.fn)
  }
  if (!is.null(sex.linker.fn)) {
    df <- add.linker.data(df, sex.linker.fn, "sex")
  }
  if (!is.null(external.id.linker.fn)) {
    df <- add.linker.data(df, external.id.linker.fn, "output")
  }
  write.table(res, output.fn, row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
}

if (exists("snakemake")) {
  run.construct.linker(
    snakemake@params[["logbook"]],
    snakemake@params[["sex_linker"]],
    snakemake@params[["external_id_linker"]],
    snakemake@output[["linker"]]
  )
}
