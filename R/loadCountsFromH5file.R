#' Gets chunk from GSE identifiers.
#' @param samples, containing a list of samples
#' @return list of chunks
#'

gsmToChunk <- function(samples) {
  chunks <- ifelse(nchar(samples) >= 7, substring(samples, 1, nchar(samples) - 4), "GSM0")
  return(paste0(chunks,"nnnn"))
}

#' @import rhdf5client
getSamples <- function(h5f, samples_id) {
  dsamples <- HSDSDataset(h5f, samples_id)
  cnt <- 1
  sampleIndexes <- list()
  while (cnt < dsamples@shape) {
    sampleIndexes <- c(sampleIndexes, dsamples[cnt:min(c(cnt + 100000, dsamples@shape))])
    cnt <- min(c(cnt + 100001, dsamples@shape))
  }
  return(unlist(sampleIndexes))
}

#' Load count matrix from remote HDF5-file
#' @param es, containing ExpressionSet loaded from GEO. Contains empty expression matrix.
#'
#' @param url, containing url of the server and root domain.
#' @param file, containing name of the file (relative to the root domain)
#' @param sampleIndexes, containing sample indexes list
#'
#' @return ExpressionSet object with loaded count matrix
#'
#' @export
#' @import data.table
#' @import rhdf5client
#'
#' @examples
#' file <- "/counts/dee2/athaliana_star_matrix_20221107.h5"
#' es <- loadCountsFromH5FileHSDS(es, url, file)
loadCountsFromH5FileHSDS <- function(es, url, file, sampleIndexes = NULL) {
  if (nrow(es) > 0) {
    return(es)
  }
  src <- httr::parse_url(url)
  dir <- src$query$domain
  src <- paste0(src$scheme, '://', src$hostname, '/', src$path)
  src <- HSDSSource(src)
  absPath <- file.path(dir, file, fsep="/")
  f <- HSDSFile(src, absPath)
  name <- basename(file)

  metafilepath <- file.path(dirname(absPath), 'meta.h5', fsep="/")
  metaf <- HSDSFile(src, metafilepath)
  metads <- HSDSDataset(metaf, '/meta')
  metatable <- metads[1:metads@shape]
  h5_meta <- metatable[file_name == name]

  if (is.null(sampleIndexes)) {
    sampleIndexes <- getSamples(f, h5_meta$sample_id)
    match_accession <- match(es$geo_accession, sampleIndexes)
    phenoData(es) <- phenoData(es[, !is.na(match_accession)])
    sampleIndexes <- match(es$geo_accession, sampleIndexes)
  }


  gene_id <- strsplit(h5_meta$gene_id, split = ":")[[1]]

  dg <- HSDSDataset(f, gene_id[[2]])
  genes <- dg[1:dg@shape]




  smap <- data.frame(sampleIndexes, geo_accession=es$geo_accession)
  smap <- smap[order(smap$sampleIndexes),]
  smap <- smap[!is.na(smap$sampleIndexes),]

  h5Indexes = list(smap$sampleIndexes,
                   seq_len(length(genes)))


  expression <- NULL
  ds <- HSDSDataset(f, '/data/expression')

  if (h5_meta$sample_dim == "rows"){
    expression <- ds[h5Indexes[[2]], h5Indexes[[1]]]
  } else {
    expression <- ds[h5Indexes[[1]], h5Indexes[[2]]]
    expression <- t(expression)
  }

  rownames(expression) <- genes
  colnames(expression) <- smap$geo_accession

  es <- es[,es$geo_accession %in% colnames(expression)]
  expression <- expression[, es$geo_accession]

  es2 <- ExpressionSet(assayData = expression,
                           phenoData = phenoData(es[, !is.na(sampleIndexes)]),
                           annotation = annotation(es),
                           experimentData = experimentData(es))
  experimentData(es2)@preprocessing$gene_counts_source <- file

  genes_annot <- strsplit(h5_meta$genes_annot, split = ";")[[1]]
  genes_annot <- unlist( lapply(strsplit(genes_annot, split = ":"), function(annot){
    setNames(annot[2], annot[1])
  }))


  genes_annot_values <- lapply(genes_annot, function(annot){
    tryCatch({
      da <- HSDSDataset(f, annot)
      as.character(da[seq_len(da@shape)])
    }, error = function(e) {})
  })
  genes_annot_values[[gene_id[1]]] <- rownames(es2)
  genes_annot_values <- genes_annot_values[!unlist(lapply(genes_annot_values, is.null))]
  fData(es2) <- cbind(fData(es2), genes_annot_values )

  return(es2)
}
#' Load count matrix from HDF5-files.
#' @param es, containing ExpressionSet loaded from GEO. Contains empty expression matrix.
#'
#' @param url, containing url of the server and root domain.
#' @return ExpressionSet with loaded count matrix
#' @export
#' @examples
#' ess <- getGEO("GSE85653")
#' es <- ess[[1]]
#' es <- loadCountsFromHSDS(es, url)
#'
loadCountsFromHSDS <- function(es, url='https://ctlab.itmo.ru/hsds/?domain=/counts') {
  if (nrow(es) > 0) {
    return(es)
  }
  src <- httr::parse_url(url)
  dir <- src$query$domain
  src <- paste0(src$scheme,'://',src$hostname,'/',src$path)
  src <- HSDSSource(src)
  priorityfilepath <- paste(dir,"/priority.h5",sep="")
  priorityf <- HSDSFile(src, priorityfilepath)
  priorityds <- HSDSDataset(priorityf, '/priority')
  priority <- data.table(priorityds[1:priorityds@shape])
  priority <- priority[, .(directory), keyby = priority]$directory

  metaindexpath <- paste(dir,"/index.h5",sep="")
  indexf <- HSDSFile(src, metaindexpath)
  sampleschunk <- unique(gsmToChunk(es$geo_accession))
  DT_counts_meta_indexes <- data.table()
  for (chunk in sampleschunk) {
      indexds <- HSDSDataset(indexf, paste0('/',chunk))
      DT_counts_meta_indexes <- rbind(DT_counts_meta_indexes, indexds[1:indexds@shape])
  }

  sample_amount <- DT_counts_meta_indexes[accession %in% es$geo_accession, .(.N), by = list(file, type_fac = factor(x = collection_type, levels = priority))]

  if (nrow(sample_amount) == 0) {
    return(es)
  }
  setorderv(x = sample_amount,cols = c("N","type_fac"),order = c(-1,1))
  destfile <- sample_amount[,.SD[1]]$file

  collection <- sample_amount[,.SD[1]]$type_fac

  DT_counts_meta_indexes <- DT_counts_meta_indexes[DT_counts_meta_indexes$file == destfile, ]

  sampleIndexes <- match(es$geo_accession, DT_counts_meta_indexes$accession)
  phenoData(es) <- phenoData(es[, !is.na(sampleIndexes)])
  sampleIndexes <- match(es$geo_accession, DT_counts_meta_indexes$accession)
  sampleIndexes <- DT_counts_meta_indexes[na.omit(sampleIndexes), ]$indexes

  es2 <- loadCountsFromH5FileHSDS(es, url, destfile, sampleIndexes)
  return(es2)
}



