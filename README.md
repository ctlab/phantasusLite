# PhantasusLite

PhantasusLite -- a tool designed to integrate the work with RNA-Seq count matrices into a single and fast R-pipeline. This R-package supports loading the data from GEO by GSE. It provides url access to the remote repository with the archs4, archs4_zoo and dee2 HDF5-files for getting the count matrix. Finally phantasusLite allows to get an ExpressionSet with the expression matrix for future differential expression analysis.

## Installation

To install the package you should:

1.  Install R-package rhdf5client from the github
2.  Install R-package phantasusLite from the github

``` r
library(devtools)
install_github("assaron/rhdf5client")
install_github("rsablina/phantasusLite")
```

## Dependencies

To run the code you need:

-   `GEOquery`
-   `rhdf5client`
-   `phantasuslite`

``` r
library(GEOquery)
library(rhdf5client)
library(phantasusLite)
```

## Quick start

To run the package enter the code sample below.

Let's load the ExpressionSet from GEO

``` r
ess <- getGEO("GSE85653")
es <- ess[[1]]
```

ExpressionSet from the GEO doesn't contain the expression matrix -- exprs(es) is empty.

``` r
head(exprs(es))
```

```         
##      GSM2280286 GSM2280287 GSM2280288 GSM2280289 GSM2280290 GSM2280291
##      GSM2280292 GSM2280293 GSM2280294 GSM2280295 GSM2280296 GSM2280297
```

Function loadCountsFromHSDS returns an ExpressionSet with the expression matrix -- the second exprs(es) contains an expression matrix.

The remote repository URL is '<https://ctlab.itmo.ru/hsds/?domain=/counts>'.

``` r
url <- 'https://ctlab.itmo.ru/hsds/?domain=/counts'
es <- loadCountsFromHSDS(es, url)
head(exprs(es))
```

```         
##           GSM2280286 GSM2280287 GSM2280288 GSM2280289 GSM2280290 GSM2280291
## AT1G01010         56         85        118         90         45         64
## AT1G01020        210        221        194        193        251        288
## AT1G01030        100         99         51         38         51         64
## AT1G01040        613        599        701       1140       1193       1415
## AT1G01050       1258       1418       1070       1624       1739       2428
## AT1G01060        270        329        291       3499       3762       5499
##           GSM2280292 GSM2280293 GSM2280294 GSM2280295 GSM2280296 GSM2280297
## AT1G01010        171        169        191         67        129        114
## AT1G01020        241        212        287        156        211        238
## AT1G01030         34         50         46         56        118         83
## AT1G01040       1188       1243       1530        506        812        837
## AT1G01050       1544       1548       1967        805       1281       1283
## AT1G01060       2585       3003       3913        225        415        355
```
