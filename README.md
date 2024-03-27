
# PhantasusLite

PhantasusLite – a tool designed to integrate the work with RNA-Seq count
matrices into a single and fast R-pipeline. This R-package supports
loading the data from GEO by GSE. It provides url access to the remote
repository with the archs4, archs4_zoo and dee2 HDF5-files for getting
the count matrix. Finally phantasusLite allows to get an ExpressionSet
with the expression matrix for future differential expression analysis.

## Installation

It is recommended to install the release version of the package from
Bioconductor using the following commands:

``` r
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("phantasusLite")
```

Alternatively, the most recent version of the package can be installed
from the GitHub repository:

``` r
library(devtools)
install_github("ctlab/phantasusLite")
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
`exprs(es)` is empty.

``` r
head(exprs(es))
```

    ##      GSM1281300 GSM1281301 GSM1281302 GSM1281303 GSM1281304 GSM1281305
    ##      GSM1281306 GSM1281307

Function loadCountsFromHSDS returns an ExpressionSet with the expression
matrix – the second exprs(es) contains an expression matrix.

The remote repository URL is
‘<https://alserglab.wustl.edu/hsds/?domain=/counts>’.

``` r
url <- 'https://alserglab.wustl.edu/hsds/?domain=/counts'
es <- loadCountsFromHSDS(es, url)
head(exprs(es))
```

    ##                    GSM1281300 GSM1281301 GSM1281302 GSM1281303 GSM1281304
    ## ENSMUSG00000000001       1015        603        561        549        425
    ## ENSMUSG00000000003          0          0          0          0          0
    ## ENSMUSG00000000028        109         34          0         14          9
    ## ENSMUSG00000000031          0         18          0          0          0
    ## ENSMUSG00000000037          0          0          0          0          0
    ## ENSMUSG00000000049          0          0          0          0          0
    ##                    GSM1281305 GSM1281306 GSM1281307
    ## ENSMUSG00000000001        853        407        479
    ## ENSMUSG00000000003          0          0          0
    ## ENSMUSG00000000028        165          0         15
    ## ENSMUSG00000000031          0          0          0
    ## ENSMUSG00000000037          0          0          0
    ## ENSMUSG00000000049          0          0          0

The available gene annotations are also filled in:

``` r
head(fData(es))
```

    ##                    Gene symbol          ENSEMBLID
    ## ENSMUSG00000000001       Gnai3 ENSMUSG00000000001
    ## ENSMUSG00000000003        Pbsn ENSMUSG00000000003
    ## ENSMUSG00000000028       Cdc45 ENSMUSG00000000028
    ## ENSMUSG00000000031         H19 ENSMUSG00000000031
    ## ENSMUSG00000000037       Scml2 ENSMUSG00000000037
    ## ENSMUSG00000000049        Apoh ENSMUSG00000000049
