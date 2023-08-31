
# PhantasusLite

PhantasusLite – a tool designed to integrate the work with RNA-Seq count
matrices into a single and fast R-pipeline. This R-package supports
loading the data from GEO by GSE. It provides url access to the remote
repository with the archs4, archs4_zoo and dee2 HDF5-files for getting
the count matrix. Finally phantasusLite allows to get an ExpressionSet
with the expression matrix for future differential expression analysis.

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

- `GEOquery`
- `rhdf5client`
- `phantasuslite`

``` r
library(GEOquery)
library(rhdf5client)
library(phantasusLite)
```

## Quick start

To run the package enter the code sample below.

Let’s load the ExpressionSet from GEO

``` r
ess <- getGEO("GSE53053")
es <- ess[[1]]
```

ExpressionSet from the GEO doesn’t contain the expression matrix –
exprs(es) is empty.

``` r
head(exprs(es))
```

    ##      GSM1281300 GSM1281301 GSM1281302 GSM1281303 GSM1281304 GSM1281305
    ##      GSM1281306 GSM1281307

Function loadCountsFromHSDS returns an ExpressionSet with the expression
matrix – the second exprs(es) contains an expression matrix.

The remote repository URL is
‘<https://ctlab.itmo.ru/hsds/?domain=/counts>’.

``` r
url <- 'https://ctlab.itmo.ru/hsds/?domain=/counts'
es <- loadCountsFromHSDS(es, url)
head(exprs(es))
```

    ##               GSM1281300 GSM1281301 GSM1281302 GSM1281303 GSM1281304 GSM1281305
    ## 0610007P14Rik         86         67         30         46         23         61
    ## 0610009B22Rik         29         22          3          0         33         13
    ## 0610009L18Rik          0          0          7          0          0         15
    ## 0610009O20Rik        103         38         17         20         31         54
    ## 0610010F05Rik        259         91        115         88        113        185
    ## 0610010K14Rik         17          6          0          0          1          0
    ##               GSM1281306 GSM1281307
    ## 0610007P14Rik        105         22
    ## 0610009B22Rik         15         26
    ## 0610009L18Rik          0          9
    ## 0610009O20Rik         24         29
    ## 0610010F05Rik        108        163
    ## 0610010K14Rik          0          7
