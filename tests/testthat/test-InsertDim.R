context("s2dv::InsertDim tests")

##############################################
  dat1 <- array(c(1:26), dim = c(dat = 1, sdate = 13, ftime = 2))
  dat2 <- array(c(1:24), dim = c(2, 3, c = 4))
##############################################
test_that("1. Input checks", {

  expect_error(
  InsertDim(c()),
  "Parameter 'data' cannot be NULL."
  )
  expect_error(
  InsertDim(list(1:10), posdim = 1, lendim = 1),
  "Parameter 'data' must be an array."
  )
  expect_error(
  InsertDim(1:10, posdim = 'a', lendim = 2),
  "Parameter 'posdim' must be a positive integer."
  )
  expect_error(
  InsertDim(1:10, posdim = 0, lendim = 2),
  "Parameter 'posdim' must be a positive integer."
  )
  expect_error(
  InsertDim(1:10, posdim = 5, lendim = 2),
  "Parameter 'posdim' cannot excess the number of dimensions of parameter 'data' plus 1"
  )
  expect_error(
  InsertDim(1:10, posdim = 1, lendim = 0.2),
  "Parameter 'lendim' must be a positive integer."
  )
  expect_error(
  InsertDim(1:10, posdim = 1, lendim = T),
  "Parameter 'lendim' must be a positive integer."
  )
  expect_warning(
  InsertDim(1:10, posdim = 1, lendim = 1), 
  "The name of new dimension is not given. Set the name as 'new'."
  )
  expect_error(
  InsertDim(1:10, posdim = 1, lendim = 1, name = 1),
  "Parameter 'name' must be a character string."
  )
  expect_error(
  InsertDim(1:10, posdim = 1, lendim = 1, ncores = 'a'),
  "Parameter 'ncores' must be a positive integer."
  )
  expect_error(
  InsertDim(1:10, posdim = 1, lendim = 1, ncores = 0),
  "Parameter 'ncores' must be a positive integer."
  )

})

##############################################
test_that("2. Output checks: dat1", {

  expect_equal(
    dim(InsertDim(dat1, posdim = 1, lendim = 2)),
    dim(array(dim = c(new = 2, dat = 1, sdate = 13, ftime = 2)))
  )
  expect_equal(
    dim(InsertDim(dat1, posdim = 3, lendim = c(d = 2))),
    c(dat = 1, sdate = 13, d = 2, ftime = 2)
  )
  expect_equal(
    as.vector(InsertDim(dat1, posdim = 1, lendim = 2)[1,,,]),
    as.vector(dat1)
  )
  expect_equal(
    as.vector(InsertDim(dat1, posdim = 1, lendim = 2)[2,,,]),
    as.vector(dat1)
  )

})

##############################################
test_that("3. Output checks: dat2", {

  expect_equal(
    dim(InsertDim(dat2, posdim = 4, lendim = 2, name = 'd')),
    c(2, 3, c = 4, d = 2)
  )

  expect_equal(
    as.vector(InsertDim(dat2, posdim = 3, lendim = 1)[,,1,]),
    as.vector(dat2)
  )

  expect_equal(
    dim(InsertDim(InsertDim(dat2, posdim = 4, lendim = 2), posdim = 1, lendim = 4)),
    dim(array(dim = c(new = 4, 2, 3, c = 4, new = 2)))
  )

})

##############################################

