downloadGPL <- function(GPL, destdir = file.path(tempdir(), "geo"), 
                        mirrorPath = "https://ftp.ncbi.nlm.nih.gov") {
  GPL <- toupper(GPL)
  stub = gsub('\\d{1,3}$','nnn',GPL,perl=TRUE)
  GPLDirPath <- '%s/platforms/%s/%s/annot'
  fullGPLDirPath <- file.path(sprintf(GPLDirPath, destdir, stub, GPL))
  
  cachedFile <- file.path(fullGPLDirPath, paste0(GPL, ".annot.gz"))
  cachedSoft <- file.path(fullGPLDirPath, paste0(GPL, ".soft"))
  cachedSoftGz <- file.path(fullGPLDirPath, paste0(GPL, ".soft.gz"))
  if (file.exists(cachedFile)) {
    if (file.size(cachedFile) == 0) {
      file.remove(cachedFile)
      message(cachedFile, ' size is 0')
    } else {
      return (cachedFile)
    }
  } else if (file.exists(cachedSoft)) {
    if (file.size(cachedSoft) == 0) {
      file.remove(cachedSoft)
      message(cachedSoft, ' size is 0')
    } else {
      return (cachedSoft)
    }
  } else if (file.exists(cachedSoftGz)) {
    if (file.size(cachedSoftGz) == 0) {
      file.remove(cachedSoftGz)
      message(cachedSoftGz, ' size is 0')
    } else {
      return (cachedSoftGz)
    }
  }
  
  annotPath <- paste0(mirrorPath,  '/geo/platforms/%s/%s/annot/%s')
  annotURL <- sprintf(annotPath,stub,GPL,paste0(GPL,'.annot.gz'))
  dir.create(fullGPLDirPath, showWarnings = FALSE, recursive = TRUE)
  targetFile <- ''
  
  for (url in annotURL){
    req <- httr::HEAD(url)
    if (httr::status_code(req) != 404) {
      # annot available
      tmp <- tempfile(pattern=paste0(GPL, ".annot.gz"), tmpdir=fullGPLDirPath)
      tryCatch({
        download.file(url, tmp)
      }, error=function (e) {
        unlink(tmp)
        stop('Could not download GPL ', GPL, e)
      })
      
      file.copy(tmp, cachedFile)
      unlink(tmp)
      targetFile <- cachedFile
    }
    if (targetFile != ''){
      break
    }
  }
  
  if (targetFile == ''){
    # need submitter
    tmp <- tempfile(pattern=paste0(GPL, ".soft"), tmpdir=fullGPLDirPath)
    apiURL <- "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi"
    submitterURL <- paste(apiURL,'?targ=self&acc=',GPL,'&form=text&view=data',sep='')
    tryCatch({
      download.file(submitterURL, tmp, headers = c("accept-encoding" = "gzip"))
    }, error=function (e) {
      unlink(tmp)
      stop('Could not download GPL ', GPL, e)
    })
    
    file.rename(tmp, cachedSoftGz)
    targetFile <- cachedSoftGz
  }
  
  if (file.size(targetFile) == 0) {
    file.remove(targetFile)
    stop('Could not download GPL ', GPL, '. Zero file')
  }
  
  con <- file(targetFile, 'r')
  first_hundred_lines <- readLines(con, n=100)
  close(con)
  found<-grep('^\\^(PLATFORM|ANNOTATION)', first_hundred_lines, ignore.case = TRUE)
  if (!length(found)) {
    file.remove(targetFile)
    stop('Could not download GPL ', GPL, '.  Possible NCBI problems')
  }
  
  return(targetFile)
}

getGPLAnnotation <- function(GPL, destdir = file.path(tempdir(), "geo")) {
  filename <- downloadGPL(GPL, destdir)
  
  txt <- data.table::fread(filename,sep="")[[1]]
  if (length(txt) == 0) {
    txt <- c("table_begin", "table_end")
  }
  ret <- GEOquery:::.parseGPLTxt(txt)
  return(ret)
}

annotateFeatureData <- function (es, destdir = file.path(tempdir(), "geo")) {
  platform <- as.character(es$platform_id)[1]
  platformParsed <- getGPLAnnotation(platform, destdir)
  
  #https://github.com/seandavi/GEOquery/blob/master/R/parseGEO.R#L569
  #############################
  vmd <- Columns(platformParsed)
  dat <- Table(platformParsed)
  ## GEO uses "TAG" instead of "ID" for SAGE GSE/GPL entries, but it is apparently
  ##     always the first column, so use dat[,1] instead of dat$ID
  ## The next line deals with the empty GSE
  tmpnames=character(0)
  if(ncol(dat)>0) {
    tmpnames=as.character(dat[,1])
  }
  ## Fixed bug caused by an ID being "NA" in GSE15197, for example
  tmpnames[is.na(tmpnames)]="NA"
  rownames(dat) <- make.unique(tmpnames)
  ## Apparently, NCBI GEO uses case-insensitive matching
  ## between platform IDs and series ID Refs ???
  dat <- dat[match(tolower(rownames(es)),tolower(rownames(dat))),]
  # Fix possibility of duplicate column names in the
  # GPL files; this is prevalent in the Annotation GPLs
  rownames(vmd) <- make.unique(colnames(Table(platformParsed)))
  colnames(dat) <- rownames(vmd)
  ##############################
  
  featureData(es) <- new('AnnotatedDataFrame',data=dat,varMetadata=vmd)
  es
}

filterFeatureAnnotations <- function(es) {
  fvarsToKeep <- c()
  if ("Gene symbol" %in% fvarLabels(es)) {
    fvarsToKeep <- c(fvarsToKeep, "Gene symbol")
  } else {
    fvarsToKeep <- c(fvarsToKeep, grep("symbol",
                                       fvarLabels(es),
                                       ignore.case = TRUE,
                                       value = TRUE))
  }
  
  if ("Gene ID" %in% fvarLabels(es)) {
    fvarsToKeep <- c(fvarsToKeep, "Gene ID")
  } else if ("ENTREZ_GENE_ID" %in%  fvarLabels(es)){
    fvarsToKeep <- c(fvarsToKeep, "ENTREZ_GENE_ID")
  }else if ("ID" %in% fvarLabels(es)) {
    fvarsToKeep <- c(fvarsToKeep, "ID")
  } else {
    fvarsToKeep <- c(fvarsToKeep, grep("entrez",
                                       fvarLabels(es),
                                       ignore.case = TRUE,
                                       value = TRUE))
  }
  
  if (length(setdiff(fvarsToKeep, "ID")) == 0) {
    fvarsToKeep <- c(fvarsToKeep,
                     setdiff(colnames(featureData(es)), fvarsToKeep))
  }
  
  
  featureData(es) <- featureData(es)[, fvarsToKeep]
  
  if (!any(sapply(fData(es),
                  function(x) identical(rownames(es), as.character(x))
  ))) {
    fData(es) <- cbind("id"=rownames(es), fData(es))
  }
  
  es
}

take <- function(x, n) {
  sapply(x, function(x) {
    x[[n]]
  })
}

filterPhenoAnnotations <- function(es) {
  phenoData(es) <- phenoData(es)[,
                                 grepl("characteristics",
                                       varLabels(es),
                                       ignore.case = TRUE) |
                                   (varLabels(es) %in% c("title",
                                                         "id",
                                                         "geo_accession"
                                   ))]
  
  chr <- varLabels(es)[grepl("characteristics",
                             varLabels(es),
                             ignore.case = TRUE)]
  
  parsePData <- function(old.phenodata) {
    old.pdata <- pData(old.phenodata)
    labels <- varLabels(old.phenodata)
    
    new.pdata <- as.data.frame(matrix(NA, nrow = nrow(old.pdata), ncol = 0))
    
    
    for (i in seq_len(ncol(old.pdata))) {
      splitted <- strsplit(as.vector(old.pdata[[i]]), ':')
      lengths <- sapply(splitted, length)
      if (any(lengths != 2 & lengths != 0)) {
        new.pdata[[labels[i]]] <- old.pdata[[i]]
      } else {
        zeros <- which(lengths == 0)
        
        splitted[zeros] <- replicate(length(zeros), list(c(NA, NA)))
        
        newnames <- unique(trimws(take(splitted, 1)))
        newnames <- newnames[which(!is.na(newnames))]
        
        
        for (j in seq_along(newnames)) {
          name <- newnames[j]
          if (!(name %in% names(new.pdata))) {
            new.pdata[[name]] <- replicate(nrow(new.pdata), NA)
          }
          indices <- which(name == trimws(take(splitted, 1)))
          new.pdata[[name]][indices] <- trimws(take(splitted, 2)[indices])
        }
      }
    }
    rownames(new.pdata) <- rownames(old.pdata)
    AnnotatedDataFrame(new.pdata)
  }
  
  if (ncol(es) > 0) {
    phenoData(es) <- parsePData(phenoData(es))
  }
  
  es
}
