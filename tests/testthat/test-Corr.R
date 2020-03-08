context("s2dv::Corr tests")

##############################################
  # dat1
  set.seed(1)
  exp1 <- array(rnorm(240), dim = c(dataset = 1, member = 2, sdate = 5,
                                   ftime = 3, lat = 2, lon = 4))

  set.seed(2)
  obs1 <- array(rnorm(120), dim = c(dataset = 1, member = 1, sdate = 5,
                                   ftime = 3, lat = 2, lon = 4))
  set.seed(2)
  na <- floor(runif(10, min = 1, max = 120))
  obs1[na] <- NA

##############################################
test_that("1. Input checks", {

  expect_error(
  Corr(c(), c()),
  "Parameter 'exp' and 'obs' cannot be NULL."
  )
  expect_error(
  Corr(c('b'), c('a')),
  "Parameter 'exp' and 'obs' must be a numeric array."
  )
  expect_error(
  Corr(c(1:10), c(2:4)),
  paste0("Parameter 'exp' and 'obs' must be at least two dimensions ",
         "containing time_dim and memb_dim.")
  )
  expect_error(
  Corr(array(1:10, dim = c(2, 5)), array(1:4, dim = c(2, 2))),
  "Parameter 'exp' and 'obs' must have dimension names."
  )
  expect_error(
  Corr(array(1:10, dim = c(a = 2, c = 5)), array(1:4, dim = c(a = 2, b = 2))),
  "Parameter 'exp' and 'obs' must have same dimension name"
  )
  expect_error(
  Corr(exp1, obs1, memb_dim = 1),
  "Parameter 'memb_dim' must be a character string."
  )
  expect_error(
  Corr(exp1, obs1, memb_dim = 'a'),
  "Parameter 'memb_dim' is not found in 'exp' or 'obs' dimension."
  )
  expect_error(
  Corr(exp1, obs1, time_dim = c('sdate', 'a')),
  "Parameter 'time_dim' must be a character string."
  )
  expect_error(
  Corr(exp1, obs1, time_dim = 'a'),
  "Parameter 'time_dim' is not found in 'exp' or 'obs' dimension."
  )
  expect_error(
  Corr(exp1, obs1, comp_dim = c('sdate', 'ftime')),
  "Parameter 'comp_dim' must be a character string."
  )
  expect_error(
  Corr(exp1, obs1, comp_dim = 'a'),
  "Parameter 'comp_dim' is not found in 'exp' or 'obs' dimension."
  )
  expect_error(
  Corr(exp1, obs1, limits = c(1,3)),
  "Paramter 'comp_dim' cannot be NULL if 'limits' is assigned."
  )
  expect_error(
  Corr(exp1, obs1, comp_dim = 'ftime', limits = c(1)),
  paste0("Parameter 'limits' must be a vector of two positive ",
                  "integers smaller than the length of paramter 'comp_dim'.")
  )
  expect_error(
  Corr(exp1, obs1, conf.lev = -1),
  "Parameter 'conf.lev' must be a numeric number between 0 and 1."
  )
  expect_error(
  Corr(exp1, obs1, method = 1),
  "Parameter 'method' must be one of 'kendall', 'spearman' or 'pearson'."
  )
  expect_error(
  Corr(exp1, obs1, conf = 1),
  "Parameter 'conf' must be one logical value."
  )
  expect_error(
  Corr(exp1, obs1, pval = 'TRUE'),
  "Parameter 'pval' must be one logical value."
  )
  expect_error(
  Corr(exp1, obs1, ncores = 1.5),
  "Parameter 'ncores' must be a positive integer."
  )
  expect_error(
  Corr(exp = array(1:10, dim = c(sdate = 1, member = 5, a = 1)),
       obs = array(1:4, dim = c(a = 1, sdate = 2, member = 2))),
  "Parameter 'exp' and 'obs' must have same length of all dimension expect 'memb_dim'."
  )
  expect_error(
  Corr(exp = array(1:10, dim = c(sdate = 2, member = 5, a = 1)),
       obs = array(1:4, dim = c(a = 1, sdate = 2, member = 2))),
  "The length of time_dim must be at least 3 to compute correlation."
  )

})

##############################################
test_that("2. Output checks: dat1", {

  expect_equal(
    dim(Corr(exp1, obs1)$corr),
    c(n_exp = 2, n_obs = 1, dataset = 1, ftime = 3, lat = 2, lon = 4)
  )
  expect_equal(
    Corr(exp1, obs1)$corr[1:6],
    c(0.11503859, -0.46959987, -0.64113021, 0.09776572, -0.32393603, 0.27565829), 
    tolerance = 0.001
  )
  expect_equal(
    length(which(is.na(Corr(exp1, obs1)$p.val))),
    2
  )
  expect_equal(
    max(Corr(exp1, obs1)$conf.lower, na.rm = T),
    0.6332941,
    tolerance = 0.001
  )
  expect_equal(
    length(which(is.na(Corr(exp1, obs1, comp_dim = 'ftime')$corr))),
    6
  )
  expect_equal(
    length(which(is.na(Corr(exp1, obs1, comp_dim = 'ftime', limits = c(2, 3))$corr))),
    2
  )
  expect_equal(
    min(Corr(exp1, obs1, conf.lev = 0.99)$conf.upper, na.rm = TRUE),
    0.2747904,
    tolerance = 0.0001
  )
  expect_equal(
    length(Corr(exp1, obs1, conf = FALSE, pval = FALSE)),
    1
  )
  expect_equal(
    length(Corr(exp1, obs1, conf = FALSE)),
    2
  )
  expect_equal(
    length(Corr(exp1, obs1, pval = FALSE)),
    3
  )
  expect_equal(
    Corr(exp1, obs1, method = 'spearman')$corr[1:6],
    c(-0.3, -0.4, -0.6, 0.3, -0.3, 0.2)
  )
  expect_equal(
    range(Corr(exp1, obs1, method = 'spearman', comp_dim = 'ftime')$p.val, na.rm = T),
    c(0.0, 0.5),
    tolerance = 0.001
  )

})

##############################################

