context("s2dv::Eno tests")

##############################################
  set.seed(1)
  dat1 <- array(rnorm(800), dim = c(dataset = 1, member = 2, sdate = 4, 
                                    ftime = 4, lat = 10, lon = 10))
  set.seed(1)
  na <- floor(runif(40, min = 1, max = 800))
  dat1[na] <- NA


  dat2 <- array(c(-5, -7, -10:10, 12, 11, 7, 16), 
                dim = c(date = 13, ftime = 2))

##############################################
test_that("1. Input checks", {

  expect_error(
  Eno(c()),
  "Parameter 'data' cannot be NULL."
  )
  expect_error(
  Eno(data = 'a'),
  "Parameter 'data' must be a numeric array."
  )
  expect_error(
  Eno(data = array(1:10, dim = c(2,5))),
  "Parameter 'data' must have dimension names."
  )
  expect_error(
  Eno(data = 1:10, time_dim = 12),
  "Parameter 'time_dim' must be a character string."
  )
  expect_error(
  Eno(data = array(1:10, dim = c(a = 2, b = 5))),
  "Parameter 'time_dim' is not found in 'data' dimension."
  )
  expect_error(
  Eno(data = array(1:10, dim = c(a = 2, sdate = 5)), na.action = na.rm),
  "Parameter 'na.action' must be a function either na.pass or na.fail."
  )
  expect_error(
  Eno(data = c(NA,1:19), na.action = na.fail),
  paste0("Calculation fails because NA is found in paratemter 'data', ",
                "which is not accepted when ",
                "parameter 'na.action' = na.fail.")
  )
  expect_error(
  Eno(data = array(1:10, dim = c(a = 2, sdate = 5)), ncores = 0.5),
  "Parameter 'ncores' must be a positive integer."
  )

})

##############################################
test_that("2. Output checks: dat1", {
  res <- Eno(dat1)

  expect_equal(
    dim(res),
    dim(array(dim = c(dataset = 1, member = 2, ftime = 4, lat = 10, lon = 10)))
  )
  expect_equal(
    length(res[which(is.na(res))]),
    1
  )
  expect_equal(
    length(res[which(res != 4)]),
    37
  )
  expect_equal(
    mean(res, na.rm = T),
    2.768103,
    tolerance = 0.0001
  )

})

##############################################

test_that("3. Output checks: dat2", {

  expect_equal(
    Eno(dat2, time_dim = 'date'),
    array(c(6.237689, 5.683186), dim = c(ftime = 2)),
    tolerance = 0.0001
  )

})
