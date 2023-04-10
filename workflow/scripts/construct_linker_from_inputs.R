if (!require(openxlsx, quietly = TRUE)) {
  install.packages("openxlsx", repos = c(
    "https://cran.case.edu/",
    "http://lib.stat.cmu.edu/R/CRAN/",
    "https://cran.yu.ac.kr/"
  ))
  require(openxlsx, quietly = TRUE)
}
library(stringr)

construct.output.stems <- function(df) {
  ## try to construct IDs of the format PMGRCID_LSID_SQID
  output.stems <- paste(df$pmgrc, df$ls, df$sq, sep = "_")
  output.stems[is.na(df$pmgrc) | is.na(df$ls) | is.na(df$sq)] <- NA
  output.stems
}

#' Search for freetext ID remappings and apply them as needed
#'
#' The idea here is that there are evidently some freetext
#' notes along the lines of
#' "SAMPLE SWAP - this is actually PMGRC-\d+-\d+-\d" that
#' only exist in freetext in rows that are labeled by
#' _unmapped_ PMGRC ID. For downstream analysis, these ID
#' mappings must be detected and applied. To make things
#' even more confusing, there doesn't appear to be any sense
#' to which columns of the logbook this information is reported in.
#'
#' The logic in this initial implementation will likely need to be
#' made more complex when different upstream freetext field types
#' are detected.
#'
#' @param df data.frame, input logbook sheet from read.xlsx
#' @param vec character vector, vector of actual PMGRC IDs from
#' current sheet
#' @return character vector, mapped version of input character vector
#' with ID relabelings applied
apply.id.mappings <- function(df, vec) {
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
  sheet.names <- openxlsx::getSheetNames(input.fn)
  pmgrc.id <- c()
  jira.tickets <- c()
  sq.id <- c()
  ls.id <- c()
  ru.id <- c()
  sex.id <- c()
  for (sheet.name in sheet.names) {
    df <- openxlsx::read.xlsx(input.fn, sheet = sheet.name, check.names = FALSE)
    colnames(df) <- tolower(colnames(df))
    if (colnames(df)[1] == "pmgrc.id") {
      ## resolve chaos
      pmgrc.col <- 1
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
      df[, pmgrc.col] <- apply.id.mappings(df, df[, pmgrc.col])
      for (row.num in seq_len(nrow(df))) {
        pmgrc.id <- c(pmgrc.id, df[row.num, pmgrc.col])
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
    pmgrc = pmgrc.id,
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
  linker.df <- read.table(linker.fn, header = TRUE, sep = "\t", quote = "", comment.char = "", stringsAsFactors = FALSE)
  rownames(linker.df) <- linker.df[, 1]
  ## handle the situation where people not yet present are being added here
  new.subjects <- linker.df[!(linker.df[, 1] %in% df$pmgrc), 1]
  if (length(new.subjects) > 0) {
    new.df <- data.frame(
      pmgrc = new.subjects,
      jira = NA,
      ru = NA,
      sq = NA,
      ls = NA,
      sex = NA,
      output = NA
    )
    df <- rbind(df, new.df)
  }
  df[, target.colname] <- linker.df[df[, "pmgrc"], 2]
  df
}

run.construct.linker <- function(logbook.fn,
                                 sex.linker.fn,
                                 external.id.linker.fn,
                                 out.fn) {
  df <- data.frame(
    pmgrc = c("A"),
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
