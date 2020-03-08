#'Add a named dimension to an array
#'
#'Insert an extra dimension into an array at position 'posdim' with length 
#''lendim'. The array repeats along the new dimension.
#'
#'@param data An array to which the additional dimension to be added.
#'@param posdim An integer indicating the position of the new dimension.
#'@param lendim An integer indicating the length of the new dimension.
#'@param name A character string indicating the name for the new dimension. 
#'  The default value is NULL.
#'@param ncores An integer indicating the number of cores to use for parallel 
#'  computation. The default value is NULL.
#'
#'@return An array as parameter 'data' but with the added named dimension.
#'
#'@examples
#'a <- array(rnorm(15), dim = c(a = 3, b = 1, c = 5, d = 1))
#'res <- InsertDim(InsertDim(a, posdim = 2, lendim = 1, name = 'e'), 4, c(f = 2))
#'dim(res)
#'
#'@import multiApply
#'@export
InsertDim <- function(data, posdim, lendim, name = NULL, ncores = NULL) {

  # Check inputs 
  ## data
  if (is.null(data)) {
    stop("Parameter 'data' cannot be NULL.")
  }
  if (is.vector(data) & !is.list(data)) {  #is vector
    data <- as.array(data)
  }
  if (!is.array(data)) {
    stop("Parameter 'data' must be an array.")
  }
  ## posdim
  if (!is.numeric(posdim)) {
    stop("Parameter 'posdim' must be a positive integer.")
  } else if (posdim %% 1 != 0 | posdim <= 0 | length(posdim) > 1) {
    stop("Parameter 'posdim' must be a positive integer.")
  }
  if (posdim > (length(dim(data)) + 1)) {
    stop("Parameter 'posdim' cannot excess the number of dimensions of parameter 'data' plus 1")
  }
  ## lendim
  if (!is.numeric(lendim)) {
    stop("Parameter 'lendim' must be a positive integer.")
  } else if (lendim %% 1 != 0 | lendim <= 0 | length(lendim) > 1) {
    stop("Parameter 'lendim' must be a positive integer.")
  }
  ## name
  if (is.null(name)) {
    if (is.null(names(lendim))) {
      name <- 'new'
      warning("The name of new dimension is not given. Set the name as 'new'.")
    } else {
      name <- names(lendim)
    }
  } else {
    if (!is.character(name) | length(name) > 1) {
      stop("Parameter 'name' must be a character string.")
    }
  }
  ## ncores
  if (!is.null(ncores)) {
    if (!is.numeric(ncores)) {
      stop("Parameter 'ncores' must be a positive integer.")
    } else if (ncores %% 1 != 0 | ncores <= 0 | length(ncores) > 1) {
      stop("Parameter 'ncores' must be a positive integer.")
    }
  }

  ###############################
  # Calculate InsertDim

  ## create output dimension
  outdim <- lendim
  if (posdim > 1) {
    outdim <- c(dim(data)[1:(posdim - 1)], outdim)
  }
  if (posdim <= length(dim(data))) {
    outdim <- c(outdim, dim(data)[posdim:length(dim(data))])
  }
  ## create output array
  outvar <- array(dim = c(outdim))
  ## give temporary names for Apply(). The name will be replaced by data in the end
  names(dim(outvar)) <- paste0('D', 1:length(outdim))
  names(dim(outvar))[posdim] <- name #'new'

  res <- Apply(data = list(outvar), 
               margins = name, #'new', 
               fun = .InsertDim, 
               dat = data,
               ncores = ncores)$output1

  if (posdim != 1) {
    if (posdim < length(outdim)) {
    res <- Reorder(res, c(1:(posdim - 1), length(outdim), posdim:(length(outdim) - 1)))
    } else {  #posdim = length(outdim)
        res <- Reorder(res, c(1:(posdim - 1), length(outdim)))
    }
  } else {
      res <- Reorder(res, c(length(outdim), 1:(length(outdim) - 1)))
  }

  return(res)
}

.InsertDim <- function(x, data) {
  x <- data
  return(x)
}
