#'Average an array along multiple dimensions
#'
#'This function returns the mean of an array along a set of dimensions and 
#'preserves the dimension names if it has.
#'
#'@param data An array to be averaged.
#'@param dims A vector of numeric or charactor string, indicating along which 
#'  dimensions to average.
#'@param na.rm A logical value indicating whether to ignore NA values (TRUE) or 
#'  not (FALSE).
#'@return An array with the same dimension as parameter 'data' except the 'dims' 
#'  dimensions. 
#'  removed.
#'
#'@examples
#'a <- array(rnorm(24), dim = c(2, 3, 4))
#'MeanDims(a, 2)
#'MeanDims(a, c(2, 3))
#'@import multiApply
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
  if (length(dims) == length(dim(data))) {
    data <- mean(data, na.rm = na.rm)
  } else {
    if (is.character(dims)) {
      dims <- which(names(dim(data)) %in% dims)
    }
    pos <- (1:length(dim(data)))[-dims]
    data <- apply(data, pos, mean, na.rm = na.rm)
  }
  return(data)
}

