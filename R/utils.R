library(rjson)
library(GEOquery)

isValidExperimentID <- function(name) {
  grepl("^((GSE|GDS)[0-9]+(-GPL[0-9]+)?)|(GPL[0-9]+)$", name, ignore.case = TRUE)
}

#'\code{getGEODir} return the directory name for the object(s) corresponding
#'     to GEO Series Identifier.
#'
#' @param name String, containing GEO identifier of the dataset.
#'     It should start with 'GSE' or 'GDS' and can include exact GPL
#'     to annotate dataset, separated with dash ('-') from the identifier.
#'
#' @param destdir Base directory.
getGEODir <- function(name, destdir = file.path(tempdir(), "geo")) {
  if (!isValidExperimentID(name)) {
    stop(name, " does not look like a valid GEO Series ID")
  }
  type <- substr(name, 1, 3)
  GEO <- unlist(strsplit(name, "-"))[1]
  stub <- gsub("\\d{1,3}$", "nnn", GEO, perl = TRUE)
  gdsDirPath <- "%s/datasets/%s/%s/soft"
  gseDirPath <- "%s/series/%s/%s/matrix"
  gplDirPath <- "%s/platforms/%s/%s/soft"
  if (type == 'GSE') {
    fullGEODirPath <- file.path(sprintf(gseDirPath, destdir, stub, GEO))
  } else if (type == "GDS") {
    fullGEODirPath <- file.path(sprintf(gdsDirPath, destdir, stub, GEO))
  } else if (type == "GPL") {
    fullGEODirPath <- file.path(sprintf(gplDirPath, destdir, stub, GEO))
  } else {
    stop("Unsupported GEO type: ", type)
  }
  dir.create(fullGEODirPath, showWarnings = FALSE, recursive = TRUE)
  
  return(fullGEODirPath)
}

# ------------------------------------------------------------------------------
# Brief data utils
# ------------------------------------------------------------------------------
parseBriefData <- function(txt) {
  tmp <- txt[grep("!\\w*?_", txt)]
  tmp <- gsub("!\\w*?_",'', tmp)
  first.eq <- regexpr(' = ', tmp)
  tmp <- cbind(substring(tmp, first = 1, last = first.eq - 1),
               substring(tmp, first = first.eq + 3))
  tmp <- tmp[tmp[,1] != "",]
  header <- split(tmp[,2],tmp[,1])
  return(header)
}


getBriefData <- function(name, destdir = tempdir()) {
  GEO <- unlist(strsplit(name, "-"))[1]
  GEOdir <- dirname(getGEODir(GEO, destdir))
  briefFile <- file.path(GEOdir, 'brief')

  if (file.exists(briefFile)) {
    message('Using cached brief file: ', briefFile)
  } else {
    url <- sprintf("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=%s&targ=self&form=text&view=brief", GEO)
    message('Trying ', url)
    resp <- httr::GET(url)
    text <- httr::content(resp, "text", "UTF-8")
    check <- grep('Could not', text)
    if (length(check) || httr::status_code(resp) != 200) {
      message('No such dataset: ', name)
      unlink(GEOdir, recursive = TRUE, force = TRUE)
      stop('Failed to download brief data on: ', GEO, '. No such dataset')
    } else {
      writeLines(text, briefFile)
      message('Stored brief data of ', GEO, ' at ', briefFile)
    }
  }
  parsedList <- parseBriefData(readLines(briefFile))
  if (length(parsedList) == 0) {
    file.remove(briefFile)
    stop('Failed to parse brief data on: ', GEO, '. Empty list')
  }
  
  return(parsedList)
}


# ------------------------------------------------------------------------------
# Utils for getting gpls data
# ------------------------------------------------------------------------------


#' Check possible annotations for GEO Dataset.
#'
#' \code{checkGPLs} returns GPL-names for
#'     the specified GEO identifier.
#'
#' @param name String, containing GEO identifier of the dataset.
#'
#' @return Vector of filenames serialized in JSON format.
#'     If there is only one GPL for that dataset, the function will
#'     return \code{name}.
#'
#' @examples
#' \dontrun{
#' checkGPLs('GSE27112')
#' checkGPLs('GSE14308')
#' }

checkGPLs <- function(name, destdir = tempdir(), mirrorPath = tempdir()) {
  spl <- unlist(strsplit(name, "-", fixed=TRUE))
  GEO <- spl[1]
  if (length(spl) == 2) {
    gpls <- jsonlite::fromJSON(checkGPLs(spl[1]))
    gpls <- intersect(gpls, name)
    return(jsonlite::toJSON(gpls))
  }
  
  cacheDir <- file.path(destdir, 'geo')
  GEOdir <- dirname(getGEODir(GEO, cacheDir))
  GPLCacheFile <- file.path(GEOdir, 'gpls')
  if (file.exists(GPLCacheFile)) {
    gpls <- readLines(GPLCacheFile)
    if (length(gpls) != 0) {
      return(jsonlite::toJSON(gpls))
    }
  }
  
  tryCatch({
    briefData <- getBriefData(name, cacheDir)
    platforms <- briefData$platform_id
    if (length(platforms) < 1) {
      stop('No platforms in brief data')
    }
    
    if (length(platforms) == 1) {
      gpls <- c(name)
    }
    
    if (length(platforms) >= 2) {
      gpls <- unlist(lapply(platforms, function (platform) {
        paste(spl[1], platform, sep='-')
      }))
    }
    writeLines(gpls, GPLCacheFile)
    return(jsonlite::toJSON(gpls))
  }, error=function(e) {
    message(paste(e, 'Trying to use checkGPLsFallback'))
    return(checkGPLsFallback(name, destdir, mirrorPath))
  })
}


checkGPLsFallback <- function(name, destdir = tempdir(), mirrorPath = tempdir()) {
  spl <- unlist(strsplit(name, "-", fixed=TRUE))
  if (length(spl) == 2) {
    gpls <- jsonlite::fromJSON(checkGPLs(spl[1]))
    gpls <- intersect(gpls, name)
    return(jsonlite::toJSON(gpls))
  }
  
  cacheDir <- file.path(destdir, 'geo')
  
  type <- substr(name, 1, 3)
  assertthat::assert_that( (type == "GDS" || type == "GSE")
                           && nchar(name) >= 4)
  
  stub <- gsub("\\d{1,3}$", "nnn", name, perl = TRUE)
  gdsurl <- "%s/geo/%s/%s/%s/"
  
  url <- sprintf(gdsurl, mirrorPath,
                 if (type == "GDS") "datasets" else "series", stub, name)
  
  cachePath <- file.path(cacheDir,
                      if (type == "GDS") "datasets" else "series",
                      stub,
                      name,
                      "gpls")
  dir.create(dirname(cachePath), recursive = TRUE, showWarnings = FALSE)
  
  if (file.exists(cachePath)) {
    gpls <- readLines(cachePath)
    if (length(gpls) != 0) {
      return(jsonlite::toJSON(gpls))
    }
  }
  
  if (type != "GDS") {
    url <- paste0(url, "matrix/")
  }
  
  gpls <- c()
  
  tryCatch({
    for (mirror in url){
      resp <- httr::GET(mirror)
      if (httr::status_code(resp) != 404) {
        break
      }
    }
    if (httr::status_code(resp) == 404) {
      warning("No such dataset")
      return(jsonlite::toJSON(c()))
    } else {
      if (type == "GDS") {
        gpls <- c(name)
      } else {
        con <- rawConnection(resp$content)
        file.names <- GEOquery:::getDirListing(con)
        close(con)
        
        file.names <- file.names[grepl(pattern = paste0("^", name),
                                       x = file.names)]
        
        file.names <- unlist(lapply(file.names, function(x) {
          paste0(substr(x, 1, regexpr("_", x) - 1))
        }))
        if (length(file.names) == 1) {
          file.names <- c(name)
        }
        gpls <- file.names
      }
      
      writeLines(gpls, cachePath)
      return(jsonlite::toJSON(gpls))
    }
  },
  error = function(e) {
    message("Problems establishing connection.")
    return(jsonlite::toJSON(c()))
  })
}


