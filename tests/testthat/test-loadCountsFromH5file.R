library(GEOquery)

test_that("loadCountsFromHSDS works correctly", {
  url <- "https://ctlab.itmo.ru/hsds/?domain=/counts"
  ess <- getGEO("GSE85653")
  es <- ess[[1]]
  es <- loadCountsFromHSDS(es, url)
  expect_true(nrow(exprs(es)) > 0)
  samples <- es$geo_accession
  es <- loadCountsFromHSDS(es, url)
  expect_true(all(samples %in% es$geo_accession))
})


test_that("loadCountsFromHSDS returns the same ExpressionSet, if it contains counts matrix", {
  url <- "https://ctlab.itmo.ru/hsds/?domain=/counts"
  ess <- getGEO("GSE10010")
  es1 <- ess[[1]]
  es2 <- loadCountsFromHSDS(es1, url)
  expect_equal(es1, es2)

})


test_that("loadCountsFromH5FileHSDS works without metadata params", {
  url <- "https://ctlab.itmo.ru/hsds/?domain=/counts"
  file <- '/counts/archs4_zoo/Arabidopsis_thaliana_count_matrix.h5'
  ess <- getGEO("GSE85653")
  es <- ess[[1]]
  es <- loadCountsFromH5FileHSDS(es, url, file)
  expect_true(nrow(exprs(es)) > 0)
})
