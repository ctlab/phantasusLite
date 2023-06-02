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
#' @param file, containing name of the file
#' @param sample_id, containing path to the dataset with sample ids
#' @param gene_id, containing path to the dataset with gene ids
#' @param gene_id_type, containing path to the dataset with gene id type
#' @param sampleIndexes, containing sample indexes list
#'
#' @return ExpressionSet object with loaded count matrix
#'
#' @export
#' @import data.table
#' @import rhdf5client
loadCountsFromH5FileHSDS <- function(es, url, file, sample_id = NULL, gene_id = NULL, gene_id_type = NULL, sampleIndexes = NULL) {
  if (nrow(es) > 0) {
    return(es)
  }
  src <- httr::parse_url(url)
  dir <- src$query$domain
  src <- paste0(src$scheme,'://',src$hostname,'/',src$path)
  src <- HSDSSource(src)
  f <- HSDSFile(src, file)
  if (is.null(sample_id) || is.null(gene_id) || is.null(gene_id_type)) {
    collection <- gsub('^.*[\\/\\]', '', dirname(file))
    metafilepath <- paste0(dir,'/', collection, '/', collection, '.h5')
    metaf <- HSDSFile(src, metafilepath)
    metads <- HSDSDataset(metaf, '/meta')
    metatable <- metads[1:metads@shape]
    if (is.null(sample_id)) {
      sample_id <- metatable$sample_id[metatable$file_name == substring(file, nchar(dir) + nchar(collection) + 3)]
    }
    if (is.null(gene_id)) {
      gene_id <- metatable$gene_id[metatable$file_name == substring(file, nchar(dir) + nchar(collection) + 3)]
    }
    if (is.null(gene_id_type)) {
      gene_id_type <- metatable$gene_id_type[metatable$file_name == substring(file, nchar(dir) + nchar(collection) + 3)]
    }
  }

  dg <- HSDSDataset(f, gene_id)
  genes <- dg[1:dg@shape]
  if (is.null(sampleIndexes)) {
    sampleIndexes <- getSamples(f, sample_id)
    match_accession <- match(es$geo_accession, sampleIndexes)
    phenoData(es) <- phenoData(es[, !is.na(match_accession)])
    sampleIndexes <- match(es$geo_accession, sampleIndexes)
  }

  datasets <- listDatasets(f)
  if ("/info/version" %in% datasets) {
    arch_version <- HSDSDataset(f, '/info/version')
  } else {
    arch_version <- HSDSDataset(f, '/meta/info/version')
  }
  arch_version <- arch_version[1:arch_version@shape[1]]
  arch_version <- as.integer(arch_version)
  if (is.na(arch_version)) {
    arch_version <- 8
  }

  ds <- HSDSDataset(f, '/data/expression')
  smap <- data.frame(sampleIndexes, geo_accession=es$geo_accession)
  smap <- smap[order(smap$sampleIndexes),]
  smap <- smap[!is.na(smap$sampleIndexes),]
  if (arch_version >= 9) {
    expression <- ds[1:ds@shape[1], smap$sampleIndexes]
  } else {
    expression <- ds[smap$sampleIndexes, 1:ds@shape[2]]
    expression <- t(expression)
  }
  rownames(expression) <- genes
  colnames(expression) <- smap$geo_accession
  es <- es[,es$geo_accession %in% colnames(expression)]
  expression <- expression[,es$geo_accession]
  es2 <- ExpressionSet(assayData = expression,
                           phenoData = phenoData(es[, !is.na(sampleIndexes)]),
                           annotation = annotation(es),
                           experimentData = experimentData(es))
  experimentData(es)@preprocessing$gene_counts_source <- file
  if (!toupper(gene_id_type) == "GENE SYMBOL") {
    tryCatch({
      gene_symbol <- HSDSDataset(f, "/meta/genes/gene_symbol")
      gene_symbol <- unlist(gene_symbol[1:gene_symbol@shape])
      fData(es2) <- cbind(fData(es2), "Gene symbol" = gene_symbol)
    }, error = function(e) {})
  }
  if (!toupper(gene_id_type) == "ENSEMBLID") {
    tryCatch({
      entrez_id <- HSDSDataset(f, "/meta/genes/ensembl_gene_id")
      entrez_id <- unlist(entrez_id[1:entrez_id@shape])
      fData(es2) <- cbind(fData(es2), "ENSEMBLID" = entrez_id)
    }, error = function(e) {})
  }

  return(es2)
}
#' Load count matrix from HDF5-files.
#' @param es, containing ExpressionSet loaded from GEO. Contains empty expression matrix.
#'
#' @param url, containing url of the server and root domain.
#'
#' @export
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

  metafilepath <- paste0(dir,'/', collection, '/', collection, '.h5')
  metaf <- HSDSFile(src, metafilepath)
  metads <- HSDSDataset(metaf, '/meta')
  metatable <- metads[1:metads@shape]
  filename <- sub('.*/', '', destfile)
  sample_id <- metatable[metatable$file_name == filename, ]$sample_id
  gene_id <- metatable[metatable$file_name == filename, ]$gene_id
  gene_id_type <- metatable[metatable$file_name == filename, ]$gene_id_type

  DT_counts_meta_indexes <- DT_counts_meta_indexes[DT_counts_meta_indexes$file == destfile, ]

  sampleIndexes <- match(es$geo_accession, DT_counts_meta_indexes$accession)
  phenoData(es) <- phenoData(es[, !is.na(sampleIndexes)])
  sampleIndexes <- match(es$geo_accession, DT_counts_meta_indexes$accession)
  sampleIndexes <- DT_counts_meta_indexes[na.omit(sampleIndexes), ]$indexes

  filename <- paste0(dir, '/', destfile)
  es2 <- loadCountsFromH5FileHSDS(es, url, filename, sample_id, gene_id, gene_id_type, sampleIndexes)
  return(es2)
}



