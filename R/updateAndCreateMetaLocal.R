#' Creates metafiles for HDF5-files
#' @param data, contains metadata
#' @param file, contains file name
#' @param dataset_name, contains dataset name
#' @return Returns NULL
#' @keywords internal
createH5 <- function(data, file, dataset_name) {
  stopifnot(requireNamespace("rhdf5"))
  if (file.exists(file)) {
    stop(sprintf("File %s already exist", file))
    # unlink(file, recursive = FALSE)
  }
  rhdf5::h5createFile(file)
  rhdf5::h5write(data, file, dataset_name)
  rhdf5::H5close()
  return(invisible(NULL))
}

#' Converts collection meta.txt files to meta.h5, putting them to the respective
#' collection folders
#' @param counts_dir, contains directory name
#' @return Returns NULL
#' @keywords internal
createMetaH5 <- function(counts_dir){
  stopifnot(requireNamespace("rhdf5"))
  collections <- list.dirs(counts_dir, full.names = FALSE)
  collections <- collections[-1]
  for (collection in collections) {
    destdir <- paste0(counts_dir, '/', collection)
    meta <- data.table()
    filename <- paste0("meta.txt")
    h5_meta <- fread(file.path(destdir, filename), index = "file_name")
    h5filename <- file.path(destdir, 'meta.h5')
    if (file.exists(h5filename)) {
      message("Skipping ", h5filename, " as it's already exists")
      next
    }
    message("Creating ", h5filename)
    createH5(h5_meta, h5filename, 'meta')
  }
  return(invisible(NULL))
}

#' Creates HDF5-File with priority
#' @param counts_dir, contains counts directory
#' @param force logical value which lets function replace existing priority file
#' @param verbose logical value which determines a content of  the output.
#' @return Returns NULL
#' @keywords internal
createPriorityH5 <- function(counts_dir, force = FALSE, verbose = FALSE){
  stopifnot(requireNamespace("rhdf5"))
  if (!dir.exists(counts_dir)) {
    message('Counts directory ', counts_dir, " does not extist" )
    return()
  }
  h5_files <- list.files(file.path(counts_dir), "\\.h5$", full.names = TRUE, recursive = TRUE)
  list_dirs <-  list.dirs(counts_dir, full.names = FALSE, recursive = TRUE)
  list_dirs <- c(".", list_dirs)
  priority_file <- file.path(counts_dir, "counts_priority.txt")
  need_create <- TRUE
  if (file.exists(priority_file)) {
    priority <- fread(priority_file)
    if (!(setequal(priority$directory,list_dirs) && length(unique(priority$priority)) == length(priority$priority))) {
      message("!!! Priority file ", priority_file , " is invalid and will be replaced")
    } else {
      need_create <- FALSE
    }
  }
  if (need_create) {
    priority <- data.table(directory = list_dirs, priority = seq_along(list_dirs))
    write.table(x = priority, file = priority_file, sep = "\t", eol = "\n", row.names = FALSE, col.names = TRUE, quote = FALSE)
  }
  createH5(priority, 'priority.h5', 'priority')
  return(invisible(NULL))
}

#' Updates indexes from HDF5-files
#' @param counts_dir, contains counts directory
#' @param force logical value which lets function replace existing index file
#' @param verbose logical value which determines a content of  the output.
#' @return Returns NULL
#' @keywords internal
updateIndexH5 <- function(counts_dir, force = FALSE, verbose = FALSE){
  stopifnot(requireNamespace("rhdf5"))
  if (!dir.exists(counts_dir)) {
    message('Counts directory ', counts_dir, " does not extist" )
    return()
  }
  meta_name <- file.path(counts_dir, "meta.rda")
  h5_files <- list.files(file.path(counts_dir), "\\.h5$", full.names = TRUE, recursive = TRUE)
  if (!length(h5_files)) {
    return()
  }
  if (!force) {
    meta_time <- as.numeric(file.mtime(meta_name))
    h5_mtime <- max(unlist(lapply(h5_files, file.mtime)))
    dirs_mtime <- lapply(file.path(counts_dir, list_dirs[-1]), file.mtime)
        if (length(dirs_mtime) > 0) {
      dir_mtime <- max(unlist(dirs_mtime))
    } else {
      dir_mtime <- -Inf
    }

    if (file.exists(meta_name) && meta_time > h5_mtime && meta_time > dir_mtime) {
      return()
    }
  }
  if (file.exists(meta_name)) {
    unlink(meta_name)
  }
  DT_counts_meta <- data.table(matrix(ncol = 3, nrow = 0, dimnames = list( NULL, c("accession", "file", "collection_type"))))
  for (cur_dir in list_dirs) {
    dir_path <- file.path(counts_dir, cur_dir)
    dir_path <- sub(pattern = "/./?$", replacement = "", x =  dir_path)
    cur_files <- list.files(path = dir_path, pattern = "\\.h5", recursive = FALSE )
    if (length(cur_files) == 0) {
      next
    }
    if (!file.exists(file.path(dir_path, "meta.txt"))) {
      if (startsWith(x = tolower(basename(dir_path)), prefix =  "archs4")) {
        updateARCHS4meta(archDir = dir_path)
      } else if (startsWith(x = tolower(basename(dir_path)), prefix = "dee2")) {
        updateDEE2meta(destDir = dir_path)
      }
    }
    message('Populating ', cur_dir , ' counts meta' )
    if (!validateCountsCollection(collectionDir = dir_path, verbose = verbose)) {
      message("!! files in ", cur_dir , " are ignored because there is not correct meta file in this directory.")
      next
    }
    DT_part <- getCountsMetaPart(counts_dir = counts_dir, collection_name = cur_dir, verbose = verbose)
    if (length(DT_part)) {
      DT_counts_meta <- rbindlist(l = list(DT_counts_meta, DT_part))
    }
    rm(DT_part)

  }
  DT_counts_meta$chunk <- gsmToChunk(DT_counts_meta$accession)
  DT_counts_meta_split <- split(DT_counts_meta, DT_counts_meta$chunk)
  createIndexH5(DT_counts_meta_split, 'index.h5')
  save(DT_counts_meta, file = meta_name, eval.promises = TRUE)
  rm(DT_counts_meta)
  return(invisible(NULL))
}


#' Gets list  with metadata
#' @param counts_dir, contains counts directory
#' @param collection_name contains name of the collection
#' @param verbose logical value which determines a content of  the output.
#' @return list with metadata
#' @keywords internal
getCountsMetaPart <- function(counts_dir, collection_name, verbose){
  destdir <- file.path(counts_dir, collection_name)
  if (!dir.exists(destdir)) {
    return()
  }
  DT_h5_meta <- data.table()
  h5_files <- list.files(destdir, "\\.h5$", full.names = FALSE)
  if (!length(h5_files)) {
    return()
  }
  h5_meta <- fread(file.path(destdir, "meta.txt"), index = "file_name")
  for (input_file in h5_files) {
    if (input_file %in% h5_meta$file_name) {
      full_name <- file.path(destdir, input_file)
      relative_path <- file.path(collection_name, input_file)
      h5f <- rhdf5::H5Fopen(full_name, flags = "H5F_ACC_RDONLY")
      accession <- rhdf5::h5read(h5f, h5_meta[file_name == input_file, ]$sample_id)
      h5_part <- data.table(accession = accession,
                            file = relative_path,
                            collection_type = collection_name, indexes = seq_along(accession))
      rhdf5::H5Fclose(h5f)
      DT_h5_meta <- rbindlist(l = list(DT_h5_meta, h5_part))
    } else {
      if (verbose) {
        message("!! ", file.path(destdir, input_file), " is ignored")
      }
    }

  }
  return(DT_h5_meta)
}

#' Writes indexes to the file
#' @param data, contains metadata
#' @param file contains the file name
#' @return Returns NULL
#' @keywords internal
createIndexH5 <- function(data, file) {
  stopifnot(requireNamespace("rhdf5"))
  rhdf5::h5createFile(file)
  names <- names(data)
  for (i in seq_along(names)) {
    rhdf5::h5write(data[[i]], file, paste0("/",names[i]))
  }
  rhdf5::h5closeAll()
  return(invisible(NULL))
}

#' Creates \code{meta.txt} file, which describes typical dee2 files.
#' @param destDir path to directory with DEE2 .h5 files.
#' @import data.table
#' @return Returns NULL
#' @keywords internal
updateDEE2meta <- function(destDir = file.path(getOption("phantasusCacheDir"), "counts/dee2")){
  dee2files <- list.files(destDir, pattern = '\\.h5$')
  DT_meta <- data.frame(matrix(ncol = 4, nrow = length(dee2files), dimnames = list(NULL, c("file_name", "sample_id", "gene_id", "gene_id_type"))))
  DT_meta$file_name <- dee2files
  DT_meta$sample_id <- "/meta/samples/geo_accession"
  DT_meta$gene_id <- "/meta/genes/ensembl_gene_id"
  genus <- vapply(strsplit(x = dee2files, split = "_"), function(x) x[1], character(1))
  for (i_file in seq_along(dee2files)) {
    if (genus[i_file] %in% c("hsapiens", "mmusculus", "drerio", "rattus")) {
      DT_meta$gene_id_type[i_file] <- "ENSEMBLID"
      next
    } else if (genus[i_file] %in% c("athaliana", "ecoli")) {
      DT_meta$gene_id_type[i_file] <- "Locus tag"

      next
    } else if (genus[i_file] %in% c("scerevisiae")) {
      DT_meta$gene_id_type[i_file] <- "Yeast id"
      next
    } else if (genus[i_file] %in% c("dmelanogaster")){
      DT_meta$gene_id_type[i_file] <- "FlyBase id"
      next
    } else if (genus[i_file] %in% c("celegans")) {
      DT_meta$gene_id_type[i_file] <- "WormBase id"
      next
    } else{
      DT_meta$gene_id_type[i_file] <- "gene"
      next
    }
  }
  write.table(x = DT_meta, file = file.path(destDir, "meta.txt"), sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
  return(invisible(NULL))
}

#' Creates \code{meta.txt} file, which describes typical archs4 and archs4Zoo files.
#' @param archDir path to directory with arch4 .h5 files.
#' @import data.table
#' @return Returns NULL
#' @keywords internal
updateARCHS4meta <- function(archDir = file.path(getOption("phantasusCacheDir"), "counts/archs4")){
  stopifnot(requireNamespace("rhdf5"))
  archs4files <- list.files(archDir, pattern = '\\.h5$')
  DT_meta <- data.frame(matrix(ncol = 4, nrow = length(archs4files), dimnames = list(NULL, c("file_name", "sample_id", 	"gene_id", "gene_id_type"))))
  DT_meta$file_name <- archs4files
  DT_meta$sample_id <- "/meta/Sample_geo_accession"
  DT_meta$gene_id <- "/meta/genes"
  genus <- tolower(vapply(strsplit(x = archs4files, split = "_"), function(x) x[1], character(1)))
  for (i_file in seq_along(archs4files)) {
    cur_file <- file.path(archDir, archs4files[i_file])
    h5f <- H5Fopen(cur_file, flags = "H5F_ACC_RDONLY")
    arch_version <- if (H5Lexists(h5f, "info/version")) {
        rhdf5::h5read(h5f, "info/version")
    } else {
        rhdf5::h5read(h5f, "meta/info/version")
    }
    arch_version <- as.integer(arch_version)
    if (is.na(arch_version)) {
      arch_version <- 8
    }
    if (arch_version >= 9) {
      DT_meta$sample_id[i_file] <- "/meta/samples/geo_accession"
      DT_meta$gene_id[i_file] <-  "/meta/genes/genes"
    }
    if (genus[i_file] %in% c("human", "mouse")) {
      DT_meta$gene_id_type[i_file] <- "Gene Symbol"
      next
    } else if (genus[i_file] %in% c("rattus", "bos", "gallus", "danio")) {
      DT_meta$gene_id_type[i_file] <- "ENSEMBLID"
      next
    } else if (genus[i_file] %in% c("arabidopsis")) {
      DT_meta$gene_id_type[i_file] <- "TAIR id"
      next
    } else if (genus[i_file] %in% c("saccharomyces")) {
      DT_meta$gene_id_type[i_file] <- "ORF id"
      next
    } else if (genus[i_file] %in% c("caenorhabditis")) {
      DT_meta$gene_id_type[i_file] <- "WormBase id"
      next
    } else if (genus[i_file] %in% c("drosophila")) {
      DT_meta$gene_id_type[i_file] <- "FlyBase id"
      next
    } else{
      DT_meta$gene_id_type[i_file] <- "gene"
    }
    rhdf5::H5Fclose(h5f)
  }
  return(invisible(NULL))
}

#' Validates counts collection
#' @param collectionDir contains directory name
#' @param verbose logical value which determines a content of  the output.
#' @return false if collection is not valid
#' @keywords internal
validateCountsCollection <- function(collectionDir, verbose=FALSE){
  stopifnot(requireNamespace("rhdf5"))
  if (!file.exists(file.path(collectionDir, "meta.txt"))) {
    if (verbose) {
      message("metafile does not exist in ",  file.path(collectionDir))
    }
    return(FALSE)
  }

  h5_meta <- fread(file.path(collectionDir, "meta.txt"), index = "file_name")
  for (input_file  in h5_meta$file_name) {
    full_path <- file.path(collectionDir, input_file)
    cur_meta <- h5_meta[file_name == input_file, ]
    if (nrow(cur_meta) > 1) {
      if (verbose) {
        message("two or more rows in meta file for ", full_path )
      }
      return(FALSE)
    }
    h5f <- rhdf5::H5Fopen(full_path, flags = "H5F_ACC_RDONLY")

    tryCatch({
      is_sample_valid <- H5Lexists(h5f, name =  cur_meta$sample_id)

      if(!is_sample_valid){
        if (verbose) {
          message("can't read sample_id in ", full_path)
        }
        return(FALSE)
      }

      gene_ids <- if (rhdf5::H5Lexists(h5f, name = cur_meta$gene_id)) {
          rhdf5::h5read(h5f, name = cur_meta$gene_id)
      } else NULL

      if (length(gene_ids) == 0) {
        if (verbose) {
          message("can't read gene_id in ", full_path)
        }
        return(FALSE)
      }
      if (length(gene_ids)  != length(unique(gene_ids))) {
        if (verbose) {
          message("Non-unique gene ids in file ", full_path)
        }
        return(FALSE)
      }
    }, finally = {
        rhdf5::H5Fclose(h5f)
    })
  }
  return(TRUE)
}

