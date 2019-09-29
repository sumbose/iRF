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
devtools::install_github("karlkumbier/iRF2.0")
```
Alternatively, the package can be installed by downloading this repository and
using the command:

```r
R CMD INSTALL iRF2.0
```

You can subsequently load the package with the usual R commands:

```r
library(iRF)
```

OSX users may need to intall gfortran to compile. This can be done with the
following commands:

```r
curl -OL http://r.research.att.com/libs/gfortran-4.8.2-darwin13.tar.bz2
sudo tar fvxz gfortran-4.8.2-darwin13.tar.bz2 -C /
```


### Workflow Overview

Here is a brief description of the algorithm implemented in this package. It assumes the default behavior and is overly simplified, but should be enough to give you an general idea of what it happening under the hood.

1. Input a numeric feature matrix `x` and a response vector `y`.
2. Iteratively train `n.iter` random forests by doing...
    1. Populate the weight vector `mtry.select.prob = rep(1, ncol(x))`, which indicating the probabilty each feature would be chosen when training the random forests.
    2. Train a random forest with `x` and `y`, and save it for later use.
    3. Update `mtry.select.prob` with the Gini importance of each feature, so that the more prediction accuracy a certain feature provides, the more likely it will be selected in the next iteration.
    4. Repeat this routine `n.iter` times.
3. Find the random forest from the iteration with highest OOB accuracy, a.k.a. `rand.forest`.
4. Run Generalized RIT on `rand.forest` by calling `gRIT`, which does...
    1. Construct `read.forest` from `rand.forest` by calling `readForest`, which does...
        1. Construct `read.forest$tree.info`, a data frame where each row corresponds to a leaf node in `rand.forest`, and each column records some metadata about that leaf. This is mostly used to construct the following two matrices.
        2. Construct `read.forest$node.feature`, a numeric sparse matrix where each row corresponds to a leaf node in `rand.forest`, and each column records the split point of (the first appearance of) all features on the path to that leaf.
        3. Construct `read.forest$node.obs`, a boolean sparse matrix where each row corresponds to an observation, and each column records if that observation falls on a certain leaf in `rand.forest`. This means `rowSums(node.obs)` should be equal to `rep(ntree, nrow(x))` where `ntree` is the number of trees in each forest.
    2. Subset `read.forest`, keeping only leaves whose prediction is `rit.param$class.id` (for classification), or is over a threshold `rit.param$class.cut` (for regression).
    3. Run [Random Intersection Tree][RIT] on `read.forest$node.feature`, with weight being the precision of each leaf times its size, i.e. the number of observations fallen into that leaf. For the RIT algorithm, each row/leaf node/decision path is considered as an observation. A total of `rit.param$ntree` RITs are grown, and the union of intersections recovered by these RITs are aggregated and stored to `ints.eval` for further inspection.
    4. Calculate importance metrics for the interactions in `ints.eval` across leaf nodes of `rand.forest`.
5. Run outer layer bootstrap stability analysis on `ints.eval` by calling `stabilityScore`, which does...
    1. Generate `n.bootstrap` bootstrap samples, a.k.a. `bs.sample`, and for each sample...
        1. Fit random forests on a sample.
        2. Extract significant interactions on the fitted  forests by calling `gRIT`.
    2. Summarize interaction importance metrics across bootstrap samples.

Iterative reweighting assigns weights proportional the predictive power of a feature. As a result, component features of a significant intersection would be given more weight, and thus tend to appear earlier in the decision path. By keeping parts of high-order intersections in the path, we essentially reduce the order of these intersections. Note, however, that iterative reweighting doesn't seem to improve the accuracy of prediction.

See [Iterative random forests to discover predictive and stable high-order interactions][iRF] and [Refining interaction search through signed iterative Random Forests][s-iRF] for a much more in-depth description, but note that this code base has evolved since their publication.
    
  [RIT]: http://arxiv.org/abs/1303.6223
  [iRF]: https://www.pnas.org/content/115/8/1943
  [s-iRF]: http://arxiv.org/abs/1810.07287



