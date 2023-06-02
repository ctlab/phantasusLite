#' Creates a data table with indexes and chunks of samples in remote HDF5-files
#'
#' @param src, contains url of the server
#' @param collections, contains names of the collections
#' @return table with samples, indexes and chunks in all HDF5-files
getIndexRemote <- function(src, collections) {
  DT_h5_meta <- data.table()
  for (collection in collections) {
    filepath <- paste0('/counts/', collection, '/')
    src <- HSDSSource(src)
    metaf <- HSDSFile(src, paste0(filepath, collection, '.h5'))
    metads <- HSDSDataset(metaf, '/meta')
    h5_meta <- metads[1:metads@shape]
    for (input_file in h5_meta$file_name) {
      full_name <- paste0(filepath, input_file)
      relative_path <- file.path(collection, input_file)
      h5f <- HSDSFile(src, full_name)
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

createIndexH5REmote <- function(src) {
  collections <- c('archs4', 'archs4_zoo', 'dee2')
  DT_h5_meta <- getIndexRemote(src, collections)
  DT_h5_meta_split <- split(DT_h5_meta, DT_h5_meta$chunk)
  createIndexH5(DT_h5_meta_split, 'index.h5')
}
