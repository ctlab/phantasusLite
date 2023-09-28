#' Returns list of all HDF5-files on HSDS-server
#' @param url, containing url of the server and root domain.
#' @param directory, containing name of the directory
#'
#' @return List of all HDF5-files on the server or all files of the collection
#'
#' @export
#' @import rhdf5client
#' @examples
#' url <- 'https://ctlab.itmo.ru/hsds/?domain=/counts'
#' getHSDSFileList(url)
#'

getHSDSFileList <- function(url='https://ctlab.itmo.ru/hsds/?domain=/counts', directory = NULL) {
    src <- httr::parse_url(url)
    dir <- src$query$domain
    src <- paste0(src$scheme,'://',src$hostname,'/',src$path)
    src <- HSDSSource(src)
    hdf5FileList  <- list()
    if (is.null(directory)) {
      directories <- listDomains(src, dir)
      directories <- directories[-grep("*\\.h5$", directories)]
      directories <- gsub(paste0(dir, '/'), '', directories)
      for (directory in directories) {
        request <- paste0(src@endpoint, "/domains?domain=",
                          dir, '/', directory)
        response <- rhdf5client:::submitRequest(request)
        domains <- response[["domains"]]
        for (domain in domains) {
          if (domain$name != paste0(dir, "/", directory, '/', directory, ".h5")) {
            hdf5FileList <- append(hdf5FileList, domain$name)
          }
        }
      }
    } else {
      request <- paste0(src@endpoint, "/domains?domain=",
                        dir, '/', directory)
      response <- rhdf5client:::submitRequest(request)
      domains <- response[["domains"]]
      for (domain in domains) {
        if (domain$name != paste0(dir, "/", directory, '/', directory, ".h5")) {
          hdf5FileList <- append(hdf5FileList, domain$name)
        }
      }
    }
    hdf5FileList <- unlist(hdf5FileList)
    return(hdf5FileList)
}
