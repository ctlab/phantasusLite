#' Creates a data table with indexes and chunks of samples in remote HDF5-files
#' @param url, contains url to the root of counts files
#' @param collections, contains names of the collections
#' @return table with samples, indexes and chunks in all HDF5-files
getIndexRemote <- function(url, collections) {
  src <- httr::parse_url(url)
  dir <- src$query$domain
  if (is.null(dir)) {
    dir <- "/"
  }
  src <- paste0(src$scheme,'://',src$hostname,'/',src$path)
  src <- HSDSSource(src)

  DT_h5_meta <- data.table()
  for (collection in collections) {
    message("Processing collection ", collection)
    collection_path <- file.path(dir, collection, fsep="/")
    metaf <- HSDSFile(src, file.path(collection_path, 'meta.h5', fsep="/"))
    metads <- HSDSDataset(metaf, '/meta')
    h5_meta <- metads[1:metads@shape]
    for (input_file in h5_meta$file_name) {
      message("Processing file ", input_file)
      full_name <- file.path(collection_path, input_file, fsep="/")
      relative_path <- file.path(collection, input_file, fsep="/")
      
      h5f <- NULL
      tryCatch({
        h5f <- HSDSFile(src, full_name)  
      }, error=function(e) {
        message("Skipping due to error:")
        message(e, "\n")
      })
      
      if (is.null(h5f)) {
        next
      }
      
      accession <- getSamples(h5f, h5_meta[h5_meta$file_name == input_file, ]$sample_id)
      h5_part <- data.table(accession = accession,
                              file = relative_path,
                              collection_type = collection, indexes = seq_along(accession))
      DT_h5_meta <- rbindlist(l = list(DT_h5_meta, h5_part))
    }

  }
  DT_h5_meta$chunk <- gsmToChunk(DT_h5_meta$accession)
  return(DT_h5_meta)
}

#' Creates HDF5-file with indexes
#' @param url, contains URL to the root of counts files
#' @param collections vector of collection names to process
#' @param destfile where to put resulting index file
#' @return Returns NULL
createIndexH5Remote <- function(url, 
                                collections=c('archs4', 'dee2'), 
                                destfile="index.h5") {
  if (file.exists(destfile)) {
    stop("File ", destfile, " alsready exists")
  }
  
  DT_h5_meta <- getIndexRemote(url, collections)
  DT_h5_meta_split <- split(DT_h5_meta, DT_h5_meta$chunk)
  createIndexH5(DT_h5_meta_split, destfile)
  return(invisible(NULL))
}
