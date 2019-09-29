#include<Rcpp.h>

using namespace Rcpp;

void writePath(LogicalVector status,
               NumericVector splitVars,
               NumericVector splitPoints,
               IntegerVector leftChildIndices,
               IntegerVector rightChildIndices,
               IntegerVector nodeVarIndices,
               NumericVector paths,
               int p,
               bool firstSplit,
               int nodeIndex,
               NumericVector currentPath,
               int *ptr2offset);

// We don't use [[Rcpp::export]] since we handcraft the RcppExports files
NumericVector ancestorPath(DataFrame treeInfo,
                          IntegerVector nodeVarIndices,
                          int p,
                          int nleaf,
                          bool firstSplit) {
    NumericVector paths(nleaf * (2*p+1));

    int offset = 0;
    NumericVector currentPath(2*p);
    writePath(treeInfo["status"],
              treeInfo["split var"], treeInfo["split point"],
              treeInfo["left daughter"], treeInfo["right daughter"],
              nodeVarIndices, paths, p, firstSplit,
              0, currentPath, &offset);
    return paths;
}

void writePath(LogicalVector status,
               NumericVector splitVars,
               NumericVector splitPoints,
               IntegerVector leftChildIndices,
               IntegerVector rightChildIndices,
               IntegerVector nodeVarIndices,
               NumericVector paths,
               int p,
               bool firstSplit,
               int nodeIndex,
               NumericVector currentPath,
               int *ptr2offset) {

    if (status[nodeIndex]) {
        const int offset = *ptr2offset;
        std::copy(currentPath.begin(), currentPath.end(),
                paths.begin() + offset);
        paths[offset + 2*p] = nodeIndex + 1;
        *ptr2offset += 2*p + 1;
        return;
    }

    int nodeVarIndex = nodeVarIndices[splitVars[nodeIndex] - 1] - 1;

    int leftChild = leftChildIndices[nodeIndex] - 1;
    NumericVector leftPath = clone(currentPath);
    if (!firstSplit || leftPath[nodeVarIndex + p] == 0) {
        leftPath[nodeVarIndex] = splitPoints[nodeIndex];
    }
    writePath(status, splitVars, splitPoints,
            leftChildIndices, rightChildIndices, nodeVarIndices,
            paths, p, firstSplit, leftChild, leftPath, ptr2offset);

    int rightChild = rightChildIndices[nodeIndex] - 1;
    NumericVector rightPath = currentPath;
    if (!firstSplit || rightPath[nodeVarIndex] == 0) {
        rightPath[nodeVarIndex + p] = splitPoints[nodeIndex];
    }
    writePath(status, splitVars, splitPoints,
            leftChildIndices, rightChildIndices, nodeVarIndices,
            paths, p, firstSplit, rightChild, rightPath, ptr2offset);
}

