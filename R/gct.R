#' Reads ExpressionSet from a GCT file.
#'
#' Only versions 1.2 and 1.3 are supported.
#'
#' @param gct Path to gct file
#'
#' @return ExpressionSet object
#'
#' @examples
#' es <- readGct(system.file("extdata/testdata/gct/test.gct", package="phantasusLite"))
#' @export
readGct <- function(gct) {
    meta <- readLines(gct, n = 3)
    version <- meta[1]
    size <- as.numeric(unlist(strsplit(meta[2], "\t")))

    if (grepl("^#1.3", version)) {
        # number of column annotations = number of additional rows
        ann.col <- size[4]

        # number of row annotations = number of additional columns
        ann.row <- size[3]
    } else if (grepl("^#1.2", version)) {
        ann.col <- 0
        ann.row <- 1
    } else {
        stop("Unsupported version of gct: use 1.2 or 1.3")
    }

    colNames <- unlist(strsplit(meta[3], "\t"))
    if (grepl("/", colNames[1])) {
        rowIdField <- sub("(.*)/(.*)", "\\1", colNames[1])
        colIdField <- sub("(.*)/(.*)", "\\2", colNames[1])
    } else {
        rowIdField <- "id"
        colIdField <- "id"
    }

    colNames[1] <- rowIdField

    t <- fread(gct,
               sep="\t", col.names = colNames,
               skip = 2 + 1 + ann.col)


    rn <- t[[1]]

    if (any(duplicated(rn))) {
        warning(sprintf("duplicated row IDs: %s; they were renamed",
                        paste0(rn[head(which(duplicated(rn)))], collapse = " ")))
        rn <- make.unique(rn)
    }

    exp.cols <- (ann.row + 2):ncol(t)
    exp <- as.matrix(t[, exp.cols, with=FALSE])
    rownames(exp) <- rn

    cn <- colnames(exp)
    if (any(duplicated(cn))) {
        warning(sprintf("duplicated row IDs: %s; they were renamed",
                        paste0(cn[head(which(duplicated(cn)))], collapse = " ")))
        cn <- make.unique(cn)
        colnames(exp) <- cn
    }
    


    fdata <- makeAnnotated(t[, seq_len(ann.row + 1), with=FALSE])
    rownames(fdata) <- rn


    # parse pData
    {
        pdata.raw <- t(fread(gct, skip = 2, nrows = ann.col + 1, header=FALSE,
                   colClasses = "character"))
        pdata <- data.frame(pdata.raw[seq_len(ncol(exp)) + 1 + ann.row,],
                            stringsAsFactors = FALSE)
        colnames(pdata) <- pdata.raw[1, ]
        colnames(pdata)[1] <- colIdField
        rownames(pdata) <- colnames(exp)
        pdata <- makeAnnotated(pdata)

        res <- ExpressionSet(exp, featureData = fdata, phenoData = pdata)
    } 

    res
}

#' Saves ExpressionSet to a GCT file (version 1.3).
#'
#' @param es ExpresionSet obeject to save
#' @param file Path to output gct file
#' @param gzip Whether to gzip apply gzip-compression for the output file#'
#' @return Result of the closing file (as in `close()` function`)
#' @examples
#' es <- readGct(system.file("extdata/testdata/gct/test.gct", package="phantasusLite"))
#' out <- tempfile(fileext = ".gct.gz")
#' writeGct(es, out, gzip=TRUE)
#' @import Biobase
#' @export
writeGct <- function(es, file, gzip=FALSE) {
    if (gzip) {
        con <- gzfile(file)
    } else {
        con <- file(file)
    }
    open(con, open="w")
    writeLines("#1.3", con)

    pd <- pData(es)
    fd <- fData(es)

    ann.col <- ncol(pData(es))
    ann.row <- ncol(fData(es))
    writeLines(sprintf("%s\t%s\t%s\t%s",
                       nrow(es), ncol(es),
                       ann.row-1, ann.col-1), con)

    if (ann.col == 0 || ann.row == 0) {
        stop("There should be at least one row and one column annotation")
    }

    idCols <- c(head(colnames(fd), 1), head(colnames(pd), 1))
    idCols <- unique(idCols)
    idCols <- paste0(idCols, collapse="/")

    writeLines(paste0(c(idCols, tail(colnames(fd), -1), pd[[1]]), collapse="\t"), con)

    ann.col.table <- t(as.matrix(pd[, tail(seq_along(pd), -1), drop=FALSE]))
    ann.col.table <- cbind(
        tail(colnames(pd), -1),
        matrix(rep(NA, (ann.row-1)*(ann.col-1)), nrow=ann.col-1),
        ann.col.table)
    ann.col.table <- unname(ann.col.table)
    write.table(ann.col.table, file=con, quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
    write.table(cbind(fData(es), exprs(es)), file=con, quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
    close(con)
}

makeAnnotated <- function(data) {
    meta <- data.frame(labelDescription = colnames(data))
    rownames(meta) <- colnames(data)

    methods::new("AnnotatedDataFrame", data = data, varMeta = meta)
    methods::new("AnnotatedDataFrame", data = data, varMeta = meta)
}
