pasteInt <- function(x) {
  # Combine interaction into single string
  x <- paste(x, collapse='_')
  return(x)
}

nameInts <- function(ints, varnames, signed=TRUE) {
  # Convert interactions indicated by indices to interactions indicated by
  # variable names. Naming convention for an interaction is:
  #   <variable1(sign)>_ <variable2(sign)>_...

  varnames <- unique(varnames)
  p <- length(varnames)
  if (signed)
    signs <- lapply(ints, function(z) ifelse(z > p, '+', '-'))
  else
    signs <- ''

  # Adjust indexing to match varnames
  ints <- lapply(ints, function(z) z %% p + p * (z == p | z == 2 * p))
  ints.name <- mapply(function(i, s) nameInt(varnames, i, s), ints, signs)
  return(ints.name)
}

nameInt <- function(varnames, idx, sgn) {
  int <- paste0(varnames[idx], sgn)
  int <- paste(sort(int), collapse='_')
  return(int)
}

int2Id <- function(int, varnames.grp, signed=FALSE, split=FALSE) {
  # Determine integer index of named variable (signed or not)
  if (!split) int <- str_split(int, '_')[[1]]

  if (signed) {
    sgn <- grep('\\+$', int)
    varnames.grp <- str_remove_all(varnames.grp, '[\\+\\-]')
    int <- str_remove_all(int, '[\\+\\-]')
  }

  varnames.grp <- unique(varnames.grp)
  id <- sapply(int, function(i) which(varnames.grp == i))
  if (signed) {
    adjust <- rep(0, length(int))
    adjust[sgn] <- length(varnames.grp)
    id <- id + adjust
  }

  return(id)
}

unsign <- function(int) {
  # Remove sign indicators from interaction strings
  return(str_replace_all(as.character(int), '[-\\+]', ''))
}

intSubsets <- function(int, split=TRUE) {
  # Generate order 1, s - 1, and s subsets of an order-s interaction

  if (!split) int <- strs_plit(as.character(int), '_')[[1]]
  if (length(int) == 1) return(int)
  sub.ord <- c(1, length(int) - 1, length(int))
  subs <- lapply(sub.ord, combn, x=int, simplify=FALSE)
  subs <- unlist(subs, recursive=FALSE)
  return(subs)
}

