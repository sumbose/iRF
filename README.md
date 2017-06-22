## iterative Random Forests (iRF)

The R package `iRF` implements iterative Random Forests, a method for
iteratively growing ensemble of weighted decision trees, and detecting
high-order feature interactions by analyzing feature usage on decision paths.
This version uses source codes from the R package `randomForest` by Andy Liaw
and Matthew Weiner and the original Fortran codes by Leo Breiman and Adele
Cutler.

To download and install the package, use `devtools`

```r
library(devtools)
devtools::install_github("sumbose/iRF")
```

You can subsequently load the package with the usual R commands:

```r
library(iRF)
```

OSX users may need to intall gfortran to compile. This can be done with the
following commands:

```r
curl -O http://r.research.att.com/libs/gfortran-4.8.2-darwin13.tar.bz2
sudo tar fvxz gfortran-4.8.2-darwin13.tar.bz2 -C /
```






