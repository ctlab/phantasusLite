test_that("read.gct works (simple version)", {
    es <- read.gct(system.file("testdata/gct/test.gct", package="phantasusLite"))
    expect_equal(dim(exprs(es)), c(5, 7))
    expect_equal(colnames(es)[1], "s1")
    expect_equal(rownames(es)[1], "1415670_at")
    expect_true(all(c("Gene symbol", "Gene ID") %in% colnames(fData(es))), setdiff(c("Gene symbol", "Gene ID"), colnames(fData(es))))
    expect_true(all(c("condition") %in% colnames(pData(es))), setdiff(c("condition"), colnames(pData(es))))
})


test_that("read.gct works with gct #1.3", {
    t <- read.gct(system.file("testdata/gct/GSE63040.gct", package="phantasusLite"))
    expect_true(nrow(t) == 10)
    expect_true(ncol(t) == 4)
})

test_that("read.gct works with gct #1.2", {
    t <- read.gct(system.file("testdata/gct/GSE63040_gct12.gct", package="phantasusLite"))
    expect_true(nrow(t) == 10)
    expect_true(ncol(t) == 4)
})

test_that("read.gct works with duplicate row names", {
    gctFile <- system.file("testdata/gct/test_dup.gct", package="phantasusLite")
    expect_warning(t <- read.gct(gctFile), "duplicate")
    expect_true(!is.null(t))
})

test_that("write.gct and read.gct work", {
    gctFile <- system.file("testdata/gct/test.gct", package="phantasusLite")
    es <- read.gct(gctFile)
    gctFile2 <- tempfile(fileext = ".gct")
    write.gct(es, gctFile2)
    es2 <- read.gct(gctFile2)

    expect_identical(exprs(es), exprs(es2))
    expect_identical(fData(es), fData(es2))
    expect_identical(pData(es), pData(es2))
})

test_that("write.gct and read.gct support id column names", {
    gctFile <- system.file("testdata/gct/test_ids.gct", package="phantasusLite")
    es <- read.gct(gctFile)
    expect_equal(colnames(fData(es))[1], "row_id")
    expect_equal(colnames(pData(es))[1], "col_id")
    gctFile2 <- tempfile(fileext = ".gct")
    write.gct(es, gctFile2)
    es2 <- read.gct(gctFile2)
    
    expect_identical(exprs(es), exprs(es2))
    expect_identical(fData(es), fData(es2))
    expect_identical(pData(es), pData(es2))
})
