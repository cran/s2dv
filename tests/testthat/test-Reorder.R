context("s2dv::Reorder tests")

##############################################
  # dat1
  dat1 <- array(c(1:30), dim = c(dat = 1, sdate = 3, ftime = 2, lon = 5))

  # dat2
  set.seed(10)
  dat2 <- array(rnorm(10), dim = c(2, 1, 5))

  # dat3
  dat3 <- array(c(1:30), dim = c(dat = 1, 3, ftime = 2, 5))

##############################################
test_that("1. Input checks", {

  expect_error(
  Reorder(c()),
  "Parameter 'data' cannot be NULL."
  )
  expect_error(
  Reorder(c(1:3)),
  "Parameter 'data' must be an array."
  )
  expect_error(
  Reorder(data = dat1, c()),
  "Parameter 'order' cannot be NULL."
  )
  expect_error(
  Reorder(data = dat1, order = list(1,2)),
  "Parameter 'order' must be a vector of numeric or character string."
  )
  expect_error(
  Reorder(data = dat1, order = TRUE),
  "Parameter 'order' must be a vector of numeric or character string."
  )
  expect_error(
  Reorder(data = dat1, order = c(-1:2)),
  "Parameter 'order' must be positive integers."
  )
  expect_error(
  Reorder(data = dat1, order = c(1:5)),
  "Parameter 'order' exceeds the dimension length of parameter 'data'."
  )
  expect_error(
  Reorder(data = dat1, order = c('dat', 'time')),
  "Parameter 'order' do not match the dimension names of parameter 'data'."
  )
  expect_error(
  Reorder(data = dat1, order = 1:3),
  paste0("The length of parameter 'order' should be the same with the ",
         "dimension length of parameter 'data'.")
  )
  expect_error(
  Reorder(data = dat1, order = 'ftime'),
  paste0("The length of parameter 'order' should be the same with the ",
         "dimension length of parameter 'data'.")
  )

})

##############################################
test_that("2. Output checks: dat1", {

  expect_equal(
    dim(Reorder(dat1, c(2,1,4,3))),
    c(sdate = 3, dat = 1, lon = 5, ftime = 2)
  )
  expect_equal(
    dim(Reorder(dat1, c('sdate', 'dat', 'lon', 'ftime'))),
    c(sdate = 3, dat = 1, lon = 5, ftime = 2)
  )
  expect_equal(
    max(Reorder(dat1, c(2, 1, 4, 3)), na.rm = TRUE),
    30
  )

})


##############################################
test_that("3. Output checks: dat2", {

  expect_equal(
    dim(Reorder(dat2, c(2, 1, 3))),
    c(1, 2, 5)
  )

})


##############################################
test_that("4. Output checks: dat3", {

  expect_equal(
    dim(Reorder(dat3, c(4, 2, 1, 3))),
    c(5, 3, dat = 1, ftime = 2)
  )

})

##############################################
