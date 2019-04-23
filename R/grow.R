#' Grow
#' @export
grow <- function(x, ...) UseMethod("grow")

grow.default <- function(x, ...)
  stop("grow has not been implemented for this class of object")

#' Grow random forest
#' @method grow randomForest
#' @export
grow.randomForest <- function(x, how.many, ...) {
  y <- update(x, ntree=how.many)
  combine(x, y)
}
