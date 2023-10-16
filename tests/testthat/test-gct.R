test_that("readGct works (simple version)", {
    es <- readGct(system.file("extdata/testdata/gct/test.gct", package="phantasusLite"))
    expect_equal(dim(exprs(es)), c(5, 7))
    expect_equal(colnames(es)[1], "s1")
    expect_equal(rownames(es)[1], "1415670_at")
    expect_true(all(c("Gene symbol", "Gene ID") %in% colnames(fData(es))), setdiff(c("Gene symbol", "Gene ID"), colnames(fData(es))))
    expect_true(all(c("condition") %in% colnames(pData(es))), setdiff(c("condition"), colnames(pData(es))))
})


test_that("readGct works with gct #1.3", {
    t <- readGct(system.file("extdata/testdata/gct/GSE63040.gct", package="phantasusLite"))
    expect_true(nrow(t) == 10)
    expect_true(ncol(t) == 4)
})

test_that("readGct works with gct #1.2", {
    t <- readGct(system.file("extdata/testdata/gct/GSE63040_gct12.gct", package="phantasusLite"))
    expect_true(nrow(t) == 10)
    expect_true(ncol(t) == 4)
})

test_that("readGct works with duplicate row names", {
    gctFile <- system.file("extdata/testdata/gct/test_dup.gct", package="phantasusLite")
    expect_warning(t <- readGct(gctFile), "duplicate")
    expect_true(!is.null(t))
})

test_that("writeGct and readGct work", {
    gctFile <- system.file("extdata/testdata/gct/test.gct", package="phantasusLite")
    es <- readGct(gctFile)
    gctFile2 <- tempfile(fileext = ".gct")
    writeGct(es, gctFile2)
    es2 <- readGct(gctFile2)

    expect_identical(exprs(es), exprs(es2))
    expect_identical(fData(es), fData(es2))
    expect_identical(pData(es), pData(es2))
})

test_that("writeGct and readGct support id column names", {
    gctFile <- system.file("extdata/testdata/gct/test_ids.gct", package="phantasusLite")
    es <- readGct(gctFile)
    expect_equal(colnames(fData(es))[1], "row_id")
    expect_equal(colnames(pData(es))[1], "col_id")
    gctFile2 <- tempfile(fileext = ".gct")
    writeGct(es, gctFile2)
    es2 <- readGct(gctFile2)

    expect_identical(exprs(es), exprs(es2))
    expect_identical(fData(es), fData(es2))
    expect_identical(pData(es), pData(es2))
})
