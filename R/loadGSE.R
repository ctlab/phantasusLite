PROTOBUF_LAYOUT_VERSION = c(0x00, 0x02)

library("rhdf5")
library("rjson")
library("GEOquery")
library("data.table")

#' Checks GSE to be supported
#' @param name GSE id, with optional GPL specification
#' @param destDir path to cache directory
#' @param combine function on how to combine results, when multiple platforms are present
#' @return logical vector if the dataset is supported or not
#' @keywords internal

checkGSEType <- function (name, destDir, combine=any) {
  spl <- unlist(strsplit(name, "-", fixed=TRUE))
  GEO <- spl[1]
  
  briefData <- getBriefData(name, destDir)
  
  gpls <- spl[2]
  if (is.na(gpls)) {
    gpls <- briefData$platform_id
  }
  gplsOK <- sapply(gpls, function(gpl) {
    gplBrief <- getBriefData(gpl, destdir=destDir)
    return(as.numeric(gplBrief$data_row_count) <= 100000)
  })
  
  return(combine(gplsOK))
}

writeToList <- function(es) {
  data <- as.matrix(exprs(es))
  colnames(data) <- NULL
  row.names(data) <- NULL
  
  pdata <- pData(es)
  row.names(pdata) <- NULL
  pdata <- as.list(pdata)
  
  rownames <- rownames(es)
  
  fdata <- fData(es)
  row.names(fdata) <- NULL
  fdata <- as.list(fdata)
  
  
  ed <- experimentData(es)
  experimentList <- as.list(expinfo(ed))
  experimentList$other <- as.list(ed@other)
  experimentList$pubMedIds <- pubMedIds(ed)
  
  res <- list(data = data, pdata = pdata, fdata = fdata,
              rownames = rownames,
              colMetaNames = varLabels(es),
              rowMetaNames = fvarLabels(es),
              experimentData = experimentList)
  res
}

#' Load GEO Dataset.
#'
#' \code{loadGEO} returns the file with serialized ExpressionSet using
#'     ProtoBuf, parsed from data downloaded from GEO by identifier.
#'
#' @param name String, containing GEO identifier of the dataset.
#'     It should start with 'GSE' or 'GDS' and can include exact GPL
#'     to annotate dataset, separated with dash ('-') from the identifier.
#'
#' @param type Type of the dataset: 'GSE' or 'GDS'. If not specified,
#'     the function will take first three letters
#'     of \code{name} variable as type.
#'
#' @return File with ProtoBuf-serialized ExpressionSet-s
#'     that were downloaded by this identifier.
#'     For GSE-datasets there can be multiple annotations, so in file will be a
#'     list mapping name with GPL to ExpressionSet.
#'
#' @examples
#' \dontrun{
#'     loadGEO("GSE27112")
#'     loadGEO("GDS4922")
#' }
#'
#' @importFrom utils tail
#' @import Biobase
#' @import GEOquery

loadGEO <- function(name, type = NA, destdir = tempdir(), 
                    mirrorPath = "https://ftp.ncbi.nlm.nih.gov", 
                    counts_path = file.path(tempdir(),"counts")) {
  cacheDir <- file.path(destdir, "geo")
  
  if (!dir.exists(cacheDir)) {
    dir.create(cacheDir)
  }
  
  geoDir <- getGEODir(name, cacheDir)
  binaryName <- paste0(name, '.bin')
  filePath <- file.path(geoDir, binaryName)
  urls <- c()
  

  ess <- getES(name, type, destdir = cacheDir,
               mirrorPath = mirrorPath, counts_path = counts_path)

  files <- list()
  assign("es", utils::tail(ess, n = 1)[[1]], envir = parent.frame())

  current_rda <- file.path(geoDir, paste0(name, ".rda"))
  bin_is_valid <- checkBinValidity(filePath, file.info(current_rda)$ctime)
  if (!bin_is_valid) {
    for (i in seq_along(ess)) {
      seriesName <- if (!grepl(pattern = "-", name) && length(ess) > 1)
        paste0(name, "-", annotation(ess[[i]])) else name
      files[[seriesName]] <- writeToList(ess[[i]])
    }
    tempBinFile <- tempfile(paste0(binaryName, ".binsave"), tmpdir = geoDir)
    protolite::serialize_pb(list(layout_version = as.raw(PROTOBUF_LAYOUT_VERSION), ess = files), tempBinFile)
    message('Saved binary file: ', tempBinFile)
    file.rename(tempBinFile, filePath)
  }


  urls <- c(urls, paste0("/geo/",unlist(strsplit(filePath, cacheDir, fixed = TRUE))[2]))
  jsonlite::toJSON(urls)
}


#' Load ExpressionSet by GEO identifier
#'
#'\code{getES} return the ExpressionSet object(s) corresponding
#'     to GEO identifier.
#'
#' @param name String, containing GEO identifier of the dataset.
#'     It should start with 'GSE' or 'GDS' and can include exact GPL
#'     to annotate dataset, separated with dash ('-') from the identifier.
#'
#' @param type Type of the dataset: 'GSE' or 'GDS'. If not specified,
#'     the function will take first three letters
#'     of \code{name} variable as type.
#'
#' @param destdir Directory for caching loaded Series and GPL
#'     files from GEO database.
#'
#' @param mirrorPath URL string which specifies the source of matrices.
#'
#' @return List of ExpressionSet objects, that were available by given
#'     in \code{name} variable GEO identifier.
#'
#' @examples
#' \dontrun{
#'     getES('GSE14308', type = 'GSE', destdir = 'cache')
#'     getES('GSE27112')
#'     getES('GDS4922')
#' }
#'
#' @export

getES <- function(name, type = NA, destdir = file.path(tempdir(),"geo"), 
                  mirrorPath = "https://ftp.ncbi.nlm.nih.gov", 
                  counts_path = file.path(tempdir(),"counts")) {
  if (is.na(type)) {
    type <- substr(name, 1, 3)
  }
  geoDir <- getGEODir(name, destdir)
  possibly.cached <- file.path(geoDir, paste0(name, ".rda"))
  if (file.exists(possibly.cached)) {
    load(possibly.cached)
    message(paste("Loaded from locally cached parsed file", possibly.cached))
  } else {
    if (type == "GSE") {
      res <- getGSE(name, destdir, mirrorPath=mirrorPath, counts_path=counts_path)
    } else {
      stop("Incorrect name or type of the dataset")
    }
    if (length(res) > 1) {
      for (i in 1:length(res)) {
        ess <- list(res[[i]])
        destfile <- file.path(geoDir,
                              paste0(name,
                                     "-",
                                     annotation(res[[i]]),
                                     ".rda"))
        message(paste("Cached dataset to ", destfile))
        save(ess, file = destfile)
      }
    }
    ess <- res
    destfile <- file.path(geoDir, paste0(name, ".rda"))
    message(paste("Cached dataset to ", destfile))
    save(ess, file = destfile)
  }
  return(ess)
}


#' Load ExpressionSet from GEO Series
#'
#'\code{getGSE} return the ExpressionSet object(s) corresponding
#'     to GEO Series Identifier.
#'
#' @param name String, containing GEO identifier of the dataset.
#'     It should start with 'GSE' or 'GDS' and can include exact GPL
#'     to annotate dataset, separated with dash ('-') from the identifier.
#'
#' @param destdir Directory for caching loaded Series and GPL
#'     files from GEO database.
#'
#' @param mirrorPath URL string which specifies the source of matrices.
#'
#' @return List of ExpressionSet objects, that were available by given
#'     in \code{name} variable GEO identifier.
#'
#' @examples
#' \dontrun{
#'     getGSE('GSE14308', destdir = 'cache')
#'     getGSE('GSE27112')
#'     getGSE('GSE53986')
#' }
#'
#' @export
#' @import rhdf5

getGSE <- function(name, destdir = tempdir(), 
                   mirrorPath = "https://ftp.ncbi.nlm.nih.gov", 
                   counts_path = file.path(tempdir(),"counts")) {
  if (!isValidExperimentID(name)) {
    stop(name, " does not look like a valid GEO Series ID")
  }
  
  if (!checkGSEType(name, destdir)) {
    stop('Currently unsupported experiment type')
  }
  
  GEO <- unlist(strsplit(name, "-"))[1]
  
  stub <- gsub("\\d{1,3}$", "nnn", GEO, perl = TRUE)
  filename <- sprintf("%s_series_matrix.txt.gz", name)
  gseurl <- "%s/geo/series/%s/%s/matrix/%s"
  
  fullGEODirPath <- getGEODir(name, destdir)
  destfile <- file.path(fullGEODirPath, filename)
  
  infile <- file.exists(destfile)
  
  # Several GEO mirrors supported, but use only one for now.
  mirrorPath <- head(mirrorPath,1)
  if (!infile) {
    tempDestFile <- tempfile(paste0(filename, ".load"), tmpdir=fullGEODirPath)
    for (url in sprintf(gseurl, mirrorPath,
                        stub, GEO, filename)){
      tryCatch({
        utils::download.file(url,
                             destfile = tempDestFile,
                             method="libcurl")
        file.rename(tempDestFile, destfile)
        infile <- TRUE
      },
      error = function(e) {
        file.remove(tempDestFile)
      },
      warning = function(w) {
        file.remove(tempDestFile)
      })
      if(file.exists(destfile)){
        message(paste("Downloaded from the", url))
        break
      }
    }
    
  } else {
    message(paste("Loading from locally found file", destfile))
  }
  
  if (infile && file.size(destfile) > 0) {
    ess <- list(suppressWarnings(getGEO(filename = destfile,
                                        destdir = fullGEODirPath,
                                        getGPL = FALSE, AnnotGPL = FALSE)))
    for (i in seq_len(length(ess))) {
      ess[[i]] <- annotateFeatureData(ess[[i]], destdir)
    }
  } else {
    gpls <- fromJSON(checkGPLs(name))
    if (length(gpls) == 0) {
      stop(paste("Dataset", name, "not found"))
    }
    if (length(gpls) == 1 && gpls == name) {
      stop(paste("Can't download dataset ", name), gpls)
      
    }
    ess <- list()
    for (i in 1:length(gpls)) {
      ess[[gpls[[i]]]] <- getGSE(gpls[[i]], destdir = destdir, 
                                 mirrorPath = mirrorPath, counts_path=counts_path)[[1]]
    }
    return(ess)
  }
  ess <- lapply(ess, filterFeatureAnnotations)
  
  useHSDS <- isHSDS(counts_path)
  
  if (is.null(useHSDS)){
    if (dir.exists(counts_path)){
      ess <- lapply(ess, loadCounts, counts_dir = counts_path)
    }
  } else{
    if (useHSDS == TRUE){
      ess <- lapply(ess, phantasusLite::loadCountsFromHSDS, url = counts_path)
    }
  }
  
  ess <- lapply(ess, filterPhenoAnnotations)
  ess <- lapply(ess, inferCondition)
  
  return(ess)
}


checkBinValidity <- function(filePath, valid_from) {
  if (!file.exists(filePath)) {
    return(FALSE)
  }
  
  if (file.info(filePath)$ctime < valid_from) {
    return(FALSE)
  }
  
  raw_proto_version <- as.raw(PROTOBUF_LAYOUT_VERSION)
  bin_ref <- protolite::serialize_pb(list(raw_proto_version))
  file_head <- readBin(con = filePath,what = raw(), n = length(bin_ref))
  
  if (!all(bin_ref == file_head)) {
    return(FALSE)
  }
  
  return(TRUE)
}