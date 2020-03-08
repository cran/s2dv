context("s2dv::MeanDims tests")

##############################################
  # dat1
  dat1 <- array(c(1:20), 
                dim = c(dat = 1, sdate = 5, ftime = 4))
  # dat2
  dat2 <- dat1
  set.seed(1)
  na <- floor(runif(4, min = 1, max = 20))
  dat2[na] <- NA

  # dat3
  set.seed(2)
  dat3 <- array(rnorm(48), dim = c(member = 2, time = 4, 2, lon = 3))

  # dat4
  set.seed(3)
  dat4 <- rnorm(10)

##############################################
test_that("1. Input checks", {

  expect_error(
  MeanDims(c()),
  "Parameter 'data' cannot be NULL."
  )
  expect_error(
  MeanDims(data = 'a'),
  "Parameter 'data' must be a numeric array."
  )
  expect_error(
  MeanDims(c(1:10), c()),
  "Parameter 'dims' cannot be NULL."
  )
  expect_error(
  MeanDims(dat1, dims = list(1,2)), 
  "Parameter 'dims' must be a vector of numeric or character string."
  )
  expect_error(
  MeanDims(dat1, dims = c(TRUE, TRUE)),
  "Parameter 'dims' must be a vector of numeric or character string."
  )
  expect_error(
  MeanDims(dat1, dims = c(0, 1)),
  "Parameter 'dims' must be positive integers."
  )
  expect_error(
  MeanDims(dat1, dims = 5),
  "Parameter 'dims' exceeds the dimension length of parameter 'data'."
  )
  expect_error(
  MeanDims(dat1, dims = 'lat'),
  "Parameter 'dims' do not match the dimension names of parameter 'data'."
  )
  expect_error(
  MeanDims(dat1, dims = 'ftime', na.rm = na.omit),
  "Parameter 'na.rm' must be one logical value."
  )

})

##############################################
test_that("2. Output checks: dat1", {

  expect_equal(
    dim(MeanDims(dat1, dims = c(1))),
    c(sdate = 5, ftime = 4)
  )
  expect_equal(
    dim(MeanDims(dat1, dims = c(1, 3))),
    c(sdate = 5)
  )
  expect_equal(
    dim(MeanDims(dat1, dims = c('sdate', 'ftime'))),
    c(dat = 1)
  )
  expect_equal(
  MeanDims(dat1, dims = c('sdate'))[1:2],
  c(3, 8)
  )

})

##############################################
test_that("3. Output checks: dat2", {

  expect_equal(
    mean(MeanDims(dat2, dims = c(1,3), na.rm = TRUE), na.rm = TRUE),
    10
  )
  expect_equal(
    mean(MeanDims(dat2, dims = c(1,3), na.rm = FALSE), na.rm = TRUE),
    11.16667,
    tolerance = 0.0001
  )
  expect_equal(
    length(which(is.na(MeanDims(dat2, dims = c(1,3), na.rm = FALSE)))),
    2
  )

})

##############################################
test_that("4. Output checks: dat3", {

  expect_equal(
    dim(MeanDims(dat3, dims = c(1))),
    c(time = 4, 2, lon = 3)
  )
  expect_equal(
    dim(MeanDims(dat3, dims = c('time', 'lon'))),
    c(member = 2, 2)
  )
  expect_equal(
    dim(MeanDims(dat3, dims = c('time', 'lon', 'member'))),
    c(a = 2, 2)[2]
  )

})

##############################################
test_that("5. Output checks: dat4", {

  expect_equal(
    length(MeanDims(dat4, dims = 1)),
    1
  )
})
##############################################

