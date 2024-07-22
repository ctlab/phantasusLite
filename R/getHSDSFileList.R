listDomainFiles <- function(src, directory) {
    request <- paste0(src@endpoint, "/domains?domain=", directory)
    response <- rhdf5client:::submitRequest(request)
    domains <- response[["domains"]]
    
    res <- lapply(domains, function(d) {
        if (d$class == "domain") {
            return(d$name)
        }
        
        stopifnot(d$class == "folder")
        return(c(d$name,
                 listDomainFiles(src, d$name)))
    })
    res <- unlist(res)    
    res
}


#' Returns list of all HDF5-files on HSDS-server
#' @param url, containing url of the server and root domain.
#' @param directory, containing name of the directory
#'
#' @return List of all HDF5-files on the server or all files of the collection
#'
#' @export
#' @import rhdf5client
#' @examples
#' url <- 'https://alserglab.wustl.edu/hsds/?domain=/counts'
#' getHSDSFileList(url)
#'

getHSDSFileList <- function(url='https://alserglab.wustl.edu/hsds/?domain=/counts', directory = NULL) {
    src <- httr::parse_url(url)
    root <- src$query$domain
    src <- paste0(src$scheme,'://',src$hostname,'/',src$path)
    src <- HSDSSource(src)
    if (is.null(directory)) {
        directory <- root
    } else {
        directory <- paste0(root, "/", directory)
    }
    
    hdf5FileList  <- listDomainFiles(src, directory)
    hdf5FileList <- grep(".*\\.h5$", hdf5FileList, value = TRUE)
    return(hdf5FileList)
}
