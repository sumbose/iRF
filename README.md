## iterative Random Forests (iRF)

The R package `iRF` implements iterative Random Forests, a method for
iteratively growing an ensemble of weighted decision trees, and detecting
high-order feature interactions by analyzing feature usage on decision paths.
This package uses source codes from the R packages `FSInteract` by Hyun Jik Kim
and Rajen D. Shah, `randomForest` by Andy Liaw and Matthew Wiener, and the
original Random Forest Fortran codes by Leo Breiman and Adele Cutler. It was
built and tsted in OSX and linux.

To download and install the package, use `devtools`

```r
install.packages('devtools')
library(devtools)
devtools::install_github("sumbose/iRF")
```

You can subsequently load the package with the usual R commands:

```r
library(iRF)
```
For a detailed description on the usage of `iRF`, see the
[vignette](https://cdn.rawgit.com/sumbose/iRF/master/vignettes/vignette2.html). 

### System requirements
This package requires C++11 and  R (>= 3.1.2). OSX users may need to install gfortran 
to compile, which can be done with the following commands:

```r
curl -O http://r.research.att.com/libs/gfortran-4.8.2-darwin13.tar.bz2
sudo tar fvxz gfortran-4.8.2-darwin13.tar.bz2 -C /
```
