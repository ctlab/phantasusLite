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

getGSE <- function(name, destdir = tempdir(), mirrorPath = "https://ftp.ncbi.nlm.nih.gov", 
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
      ess[[gpls[[i]]]] <- getGSE(gpls[[i]], destdir = destdir, mirrorPath = mirrorPath)[[1]]
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

