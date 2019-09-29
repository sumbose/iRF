#' getTree
#'
#' Read out tree metadata from randomForest or ranger objects
#'
#' @param rfobj an object of class randomForest or ranger
#' @param k tree index to read
#' @param labelVar
#' @param terminal node membership, required for reading ranger
#'
#' @import dplyr
#' @importFrom fastmatch "%fin%"
getTree <- function(x, ...) UseMethod("getTree")

getTree.default <- function(...)
    stop(deparse(substitute(rfobj)), "is not of class ranger or randomForest")

getTree.ranger <- function(rfobj, k=1) {
  # Read metadata from forest
  tree.info <- ranger::treeInfo(rfobj, k) %>%
      transmute(`left daughter` = leftChild+1L,
                `right daughter` = rightChild+1L,
                `split var` = splitvarID+1L,
                `split point` = splitval,
                status = terminal,
                prediction = prediction)

  return(tree.info)
}

getTree.randomForest <- function(rfobj, k=1) {
  # Check whether current tree can be read
  if (is.null(rfobj$forest)) {
    stop("No forest component in ", deparse(substitute(rfobj)))
  }
  if (k > rfobj$ntree) {
    stop("There are fewer than ", k, "trees in the forest")
  }

  # Read metadata from forest
  if (rfobj$type == "regression") {
      tree.info <- data.frame(rfobj$forest$leftDaughter[,k],
                              rfobj$forest$rightDaughter[,k],
                              rfobj$forest$bestvar[,k],
                              rfobj$forest$xbestsplit[,k],
                              rfobj$forest$nodestatus[,k] == -1,
                              rfobj$forest$nodepred[,k])
  } else {
      tree.info <- data.frame(rfobj$forest$treemap[,,k],
                              rfobj$forest$bestvar[,k],
                              rfobj$forest$xbestsplit[,k],
                              rfobj$forest$nodestatus[,k] == -1,
                              rfobj$forest$nodepred[,k])
  }
  tree.info <- tree.info[1:rfobj$forest$ndbigtree[k], ]

  colnames(tree.info) <- c("left daughter", "right daughter", "split var", 
                           "split point", "status", "prediction")

  return(tree.info)
}

