## iterative Random Forests (iRF)

The R package `iRF` implements iterative Random Forests, a method for
iteratively growing ensemble of weighted decision trees, and detecting
high-order feature interactions by analyzing feature usage on decision paths.
This version uses source codes from the R package `randomForest` by Andy Liaw
and Matthew Weiner and the original Fortran codes by Leo Breiman and Adele
Cutler.

# Installation
To install `iRF` from CRAN:
```r
install.packages('iRF')
```

To install the development version from GitHub:

```r
install.packages('devtools')
devtools::install_github("sumbose/iRF")
```

You can subsequently load the package with the usual R commands:
```r
library(iRF)
```

This package requires gfortran to compile. OSX users may need to intall gfortran
available [here](https://cran.r-project.org/bin/macosx/tools/). For details on
usage, see our
[vignette](https://www.stat.berkeley.edu/~kkumbier/vignette.html).
