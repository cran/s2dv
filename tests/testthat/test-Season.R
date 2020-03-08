context("s2dv::Season tests")

##############################################
  # dat1
  set.seed(1)
  dat1 <- array(rnorm(144*3), dim = c(member = 2, sdate = 12*3, ftime = 2, lon = 3))

  # dat2
  dat2 <- dat1
  set.seed(2)
  na <- floor(runif(30, min = 1, max = 144*3))
  dat2[na] <- NA

##############################################
test_that("1. Input checks", {

  expect_error(
  Season(c()),
  "Parameter 'data' cannot be NULL."
  )
  expect_error(
    Season(c(NA, NA)),
    "Parameter 'data' must be a numeric array."
  )
  expect_error(
    Season(list(a = array(rnorm(50), dim = c(dat = 5, sdate = 10)), b = c(1:4))),
    "Parameter 'data' must be a numeric array."
  )
  expect_error(
    Season(array(1:10, dim = c(2, 5))),
    "Parameter 'data' must have dimension names."
  )
  expect_error(
    Season(dat1, time_dim = 2),
    "Parameter 'time_dim' must be a character string."
  )
  expect_error(
    Season(dat1, time_dim = 'a'),
    "Parameter 'time_dim' is not found in 'data' dimension."
  )
  expect_error(
    Season(dat1, monini = 0, moninf = 1, monsup = 2),
    "Parameter 'monini' must be a positive integer between 1 and 12."
  )
  expect_error(
    Season(dat1, monini = 1, moninf = 'jan', monsup = 2),
    "Parameter 'moninf' must be a positive integer between 1 and 12."
  )
  expect_error(
    Season(dat1, monini = 1, moninf = 1, monsup = 'Jan'),
    "Parameter 'monsup' must be a positive integer between 1 and 12."
  )
  expect_error(
    Season(dat1, monini = 1, moninf = 1, monsup = 2, method = 'mean'),
    "Parameter 'method' should be an existing R function, e.g., mean or sum."
  )
  expect_error(
    Season(dat1, monini = 1, moninf = 1, monsup = 2, na.rm = na.omit),
  "Parameter 'na.rm' must be one logical value."
  )
  expect_error(
    Season(dat1, monini = 1, moninf = 1, monsup = 2, ncores = T),
    "Parameter 'ncores' must be a positive integer."
  )

})

##############################################
test_that("2. Output checks: dat1", {

  expect_equal(
    dim(Season(dat1, monini = 1, moninf = 1, monsup = 2)),
    c(sdate = 3, member = 2, ftime = 2, lon = 3)
  )
  expect_equal(
    dim(Season(dat1, time_dim = 'lon', monini = 1, moninf = 1, monsup = 2)),
    c(lon = 1, member = 2, sdate = 36, ftime = 2)
  )
  expect_equal(
    dim(Season(dat1, monini = 10, moninf = 12, monsup = 2)),
    c(sdate = 3, member = 2, ftime = 2, lon = 3)
  )
  expect_equal(
    median(Season(dat1, monini = 10, moninf = 12, monsup = 2)),
    0.007925,
    tolerance = 0.0001
  )
  expect_equal(
    median(Season(dat1, monini = 10, moninf = 2, monsup = 5, method = sum)),
    0.2732015,
    tolerance = 0.0001
  )

})

##############################################
test_that("3. Output checks: dat2", {

  expect_equal(
    median(Season(dat2, monini = 10, moninf = 12, monsup = 2)),
    -0.01986671,
    tolerance = 0.0001
  )
  expect_equal(
    median(Season(dat2, monini = 10, moninf = 12, monsup = 2, na.rm = F), na.rm = TRUE),
    0.06207006,
    tolerance = 0.0001
  )
  res <- Season(dat2, monini = 10, moninf = 12, monsup = 2, na.rm = F)
  expect_equal(
    length(res[which(is.na(as.vector(res)))]),
    10
  )
  res <- Season(dat2, monini = 10, moninf = 12, monsup = 2)
  expect_equal(
    length(res[which(is.na(as.vector(res)))]),
    0
  )
})

##############################################

