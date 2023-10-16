#' Removes repeated words from conditions
#' @param titles, contains titles
#' @return titles without repeated words
#' @keywords internal
removeRepeatWords <- function(titles) {
  titles_without_repeat_words <- titles
  repeat_words <- regmatches(titles, regexpr("(?![-+}\\]\\)])(\\W|_)*\\w*$", titles, ignore.case = TRUE, perl = TRUE))
  lsuff <- lcSuffix(repeat_words, ignore.case = TRUE)
  if ((all(stringr::str_length(lsuff) == stringr::str_length(repeat_words))) & (all(sub("(\\W*|_)*\\w*$", "", titles, ignore.case = TRUE, perl = TRUE) != ""))) {
    titles_without_repeat_words <- sub("(?![-+}\\]\\)])(\\W|_)*\\w*$", "", titles, ignore.case = TRUE, perl = TRUE)
  }
  if (all(grepl("-$", titles_without_repeat_words, ignore.case = TRUE))) titles_without_repeat_words <- sub("-$", "", titles_without_repeat_words, ignore.case = TRUE, perl = TRUE)
  return(titles_without_repeat_words)
}

#' Creates condition from the samples titles
#' @param gse_titles, contains titles
#' @return List of conditions and replicates
#' @keywords internal
inferConditionImpl <- function(gse_titles) {
  inferCondition <- gse_titles
  rep_num <- NULL
  if ((length(inferCondition) > 40) | (length(inferCondition) < 3))
  {
    return(list())
  } else if (! all(grepl("[A-z]", inferCondition, ignore.case = TRUE)))
  {
    return(list())
  } else if (all(grepl(" vs[. ]{1}", inferCondition)))
  {
    return(list())
  } else
  {
    lsuff <- lcSuffix(inferCondition)
    suff_length <- stringr::str_length(lsuff)
    if (suff_length > 1) inferCondition <- stringr::str_sub(inferCondition, 1, -suff_length-1)
    if (all(grepl("((bio[A-z]*)|(tech[A-z]*))?[ .#_-]*((rep.*)|(set.*)|(sample)|(exp.*)|(case))[ .#_-]*\\d+", inferCondition, ignore.case = TRUE)))
    {
      sample <- regmatches(inferCondition, regexpr("((bio[A-z]*)|(tech[A-z]*))?[ .#_-]*((rep.*)|(set.*)|(sample)|(exp.*)|(case))[ .#_-]*\\d+", inferCondition, ignore.case = TRUE))
      rep_num <- regmatches(sample, regexpr("\\d+$", sample))
      inferCondition <- sub("(?![-+}\\]\\)])(\\W|_)*((bio[A-z]*)|(tech[A-z]*))?[ .#_-]*((rep.*)|(set.*)|(sample)|(exp.*)|(case))[ .#_-]*\\d+(\\W|_)*", "", inferCondition, ignore.case = TRUE, perl = TRUE)

    }
    else if (all(grepl("((bio[A-z]*)|(tech[A-z]*))?[ .#_-]*((rep.*)|(set.*)|(sample)|(exp.*)|(case))[ .#_-]*[A-z]$", inferCondition, ignore.case = TRUE)))
    {
      sample <- regmatches(inferCondition, regexpr("((bio[A-z]*)|(tech[A-z]*))?[ .#_-]*((rep.*)|(set.*)|(sample)|(exp.*)|(case))[ .#_-]*[A-z]$", inferCondition, ignore.case = TRUE))
      rep_num <- regmatches(sample, regexpr("[A-z]$", sample))
      inferCondition <- sub("(?![-+}\\]\\)])(\\W|_)*((bio[A-z]*)|(tech[A-z]*))?[ .#_-]*((rep.*)|(set.*)|(sample)|(exp.*)|(case))[ .#_-]*[A-z]$", "", inferCondition, ignore.case = TRUE, perl = TRUE)

    }
    else if (length(unique(sub("[- \\.#_]*\\d+$", "", inferCondition))) == 1) {
      return(list())
    }
    else if (all(grepl("[- \\.#_]*\\d+$", inferCondition)) & (length(unique(regmatches(inferCondition, regexpr("\\d+$", inferCondition)))) > 1))
    {
      sample <- sub("\\d+$", "", inferCondition)
      if (length(unique(regmatches(sample, regexpr("[- \\.#_]$", sample)))) > 1)
      {
        return(list())
      }
      else
      {
        rep_num <- regmatches(inferCondition, regexpr("\\d+$", inferCondition))
        inferCondition <- removeRepeatWords(sub("[ \\.#_]*\\d+$", "", inferCondition))

      }
    }
    else return(list())
  }

  if ((length(rep_num) > 1) & (length(unique(inferCondition)) < length(unique(gse_titles))) & (length(unique(inferCondition)) > 1))
    return(list(condition=inferCondition, replicate=rep_num))
  else return(list())
}

#' Adds condition to the annotation.
#' @param es, contains ExpressionSet object
#'
#' @export
#' @return Annotated ExpressionSet with conditions and replicates
#' @examples
#' ess <- GEOquery::getGEO("GSE143903")
#' es <- ess[[1]]
#' es <- inferCondition(es)
#' es$condition # contains inferred groups
#' es$replicate # contains inferred replicate numbers
#'
inferCondition <- function(es) {
  newAnnot <- inferConditionImpl(es$title)
  if (length(newAnnot) == 2) {
    pData(es)$condition <- newAnnot$condition
    pData(es)$replicate <- newAnnot$replicate
  }
  return(es)
}

