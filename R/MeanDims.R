#'Average an array along multiple dimensions
#'
#'This function returns the mean of an array along a set of dimensions and 
#'preserves the dimension names if it has.
#'
#'@details It is recommended to use \code{'apply(x, dim, mean)'} to improve the
#'  efficiency when the dimension to be averaged is only one. 
#'
#'@param data An array to be averaged.
#'@param dims A vector of numeric or charactor string, indicating along which 
#'  dimensions to average.
#'@param na.rm A logical value indicating whether to ignore NA values (TRUE) or 
#'  not (FALSE). The default value is FALSE.
#'
#'@return An array with the same dimension as parameter 'data' except the 'dims' 
#'  dimensions. 
#'  removed.
#'
#'@examples
#'a <- array(rnorm(24), dim = c(a = 2, b= 3, c = 4))
#'print(dim(MeanDims(a, 2)))
#'print(dim(MeanDims(a, c(2, 3))))
#'print(dim(MeanDims(a, c('a', 'b'))))
#'@export
MeanDims <- function(data, dims, na.rm = FALSE) {

  # Check inputs 
  ## data
  if (is.null(data)) {
    stop("Parameter 'data' cannot be NULL.")
  }
  if (!is.numeric(data)) {
    stop("Parameter 'data' must be a numeric array.")
  }
  if (is.null(dim(data))) {  #is vector, turn into array
    data <- as.array(data)
  }
  ## dims
  if (is.null(dims)) {
    stop("Parameter 'dims' cannot be NULL.")
  }
  if (!is.vector(dims) | (is.vector(dims) & !is.numeric(dims) & !is.character(dims))) {
    stop("Parameter 'dims' must be a vector of numeric or character string.")
  }
  if (is.numeric(dims)) {
    if (any(dims < 1) | any(dims %% 1 != 0)) {
      stop("Parameter 'dims' must be positive integers.")
    } else if (any(dims > length(dim(data)))) {
      stop("Parameter 'dims' exceeds the dimension length of parameter 'data'.")
    } 
  }
  if (is.character(dims)) {
    if (!all(dims %in% names(dim(data)))) {
      stop("Parameter 'dims' do not match the dimension names of parameter 'data'.")
    }
  }
  ## na.rm
  if (!is.logical(na.rm) | length(na.rm) > 1) {
    stop("Parameter 'na.rm' must be one logical value.")
  }



  ###############################
  # Calculate MeanDims

  ## Change character dims into indices
  if (is.character(dims)) {
    tmp <- rep(0, length(dims))
    for (i in 1:length(dims)) {
      tmp[i] <- which(names(dim(data)) == dims[i])
    }
  dims <- tmp
  }

  if (length(dim(data)) == 1) {
    res <- mean(data, na.rm = na.rm) 
  } else {

    margins <- setdiff(c(1:length(dim(data))), dims)
    res <- as.array(apply(data, margins, mean, na.rm = na.rm))
    if (!is.null(names(dim(data))[margins]) & is.null(names(dim(res)))) {
      names(dim(res)) <- names(dim(data))[margins]
    }
  }

  return(res)

}

